suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(ggtext)
    library(patchwork)
})

tool_cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F",
    "NanoTS"       = "#E05A8A"
)

cat_cols <- c("Other FP" = "grey70", "RNA editing (REDIportal)" = "#2166AC")

add_tech <- \(dt) {
    dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
    dt[tech == "Baylor-IsoSeq", tech := "IsoSeq"]
    dt[, cell_line := sub("^([^-]+)-.*$", "\\1", sample)]
    dt[cell_line == "HG002a", cell_line := "HG002"]
    dt <- dt[cell_line %in% c("HG004", "HG005")]
    dt[, cell_line := factor(cell_line, levels = c("HG004", "HG005"))]
    dt[, tech := ifelse(grepl("MasSeq|IsoSeq", tech),
                         paste0("<span style='color:#54278f;'>", tech, "</span>"),
                         paste0("<span style='color:#d95f0e;'>", tech, "</span>"))]
    dt
}

common_theme <- list(
    theme_classic(),
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        strip.text = element_markdown(size = 15),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 14),
        aspect.ratio = 1
    )
)

# a: A>G/T>C FPs at REDIportal sites; b/c: change in precision/recall after removing them
redi_files <- args[[1]][grepl("/rna_editing/", args[[1]]) & !grepl("/spectrum/", args[[1]])]
redi <- rbindlist(lapply(redi_files, fread))
redi <- redi[Type == "SNP" & tool != "GATK"]
redi[, coverage := as.integer(coverage)]
redi <- add_tech(redi)

spec_files <- args[[1]][grepl("/rna_editing/spectrum/", args[[1]])]
spec <- rbindlist(lapply(spec_files, fread))
spec <- spec[tool != "GATK" & !grepl(",", substitution)]
spec <- add_tech(spec)

REF_COV <- 5
fp_atogttoc <- spec[BD == "FP" & substitution %in% c("A>G", "T>C") & coverage == REF_COV,
                     .(FP_AtoG_TtoC = sum(count)), by = .(tool, cell_line, tech)]
redi_agg <- redi[coverage == REF_COV, .(FP_redi_strict = sum(FP_in_REDIportal_strict)),
                  by = .(tool, cell_line, tech)]

da <- merge(fp_atogttoc, redi_agg, by = c("tool", "cell_line", "tech"))
da[, FP_other := FP_AtoG_TtoC - FP_redi_strict]
da[, pct := round(100 * FP_redi_strict / FP_AtoG_TtoC, 1)]

long_a <- melt(da, id.vars = c("tool", "cell_line", "tech", "pct", "FP_AtoG_TtoC"),
               measure.vars = c("FP_redi_strict", "FP_other"),
               variable.name = "category", value.name = "count")
long_a[, category := factor(category, levels = c("FP_other", "FP_redi_strict"),
                             labels = c("Other FP", "RNA editing (REDIportal)"))]

pa <- ggplot(long_a, aes(tool, count, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(data = da, aes(x = tool, y = FP_AtoG_TtoC, label = paste0(pct, "%")),
              inherit.aes = FALSE, angle = 45, hjust = 0, vjust = 0, size = 3,
              position = position_nudge(x = -0.2)) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    scale_fill_manual(values = cat_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.32))) +
    common_theme +
    labs(x = "Variant Caller", y = "Number of A>G/T>C FP SNVs", fill = "Category")

hp <- fread(args[[1]][grepl("combined_SNP\\.csv$", args[[1]])])[tool != "GATK"]
hp <- add_tech(hp)
.agg <- \(x) x[, .(TP = sum(TRUTH.TP), FP = sum(QUERY.FP), FN = sum(TRUTH.FN)),
               by = .(tool, cell_line, tech, coverage)]
d <- merge(.agg(hp[set == "orig"]), .agg(hp[set == "masked"]),
           by = c("tool", "cell_line", "tech", "coverage"), suffixes = c("_orig", "_masked"))
d[, `:=`(dPrecision = TP_masked / (TP_masked + FP_masked) - TP_orig / (TP_orig + FP_orig),
         dRecall    = TP_masked / (TP_masked + FN_masked) - TP_orig / (TP_orig + FN_orig),
         dF1        = 2 * TP_masked / (2 * TP_masked + FP_masked + FN_masked) -
                      2 * TP_orig / (2 * TP_orig + FP_orig + FN_orig))]
d[, coverage := factor(coverage, levels = sort(unique(coverage)))]

.p <- \(y, lab) ggplot(d, aes(coverage, .data[[y]], color = tool, group = tool)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    scale_color_manual(values = tool_cols) +
    common_theme +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) +
    labs(x = "Coverage Cutoff (DP >= n)", y = lab, color = "Variant Caller")

gg <- pa / .p("dPrecision", expression(Delta * " Precision (REDIportal-filtered)")) /
      .p("dRecall", expression(Delta * " Recall (REDIportal-filtered)")) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 20, face = "bold"),
          legend.position = "bottom", legend.box = "vertical")

ggsave(args[[2]], gg, width = 28, height = 40, units = "cm")
