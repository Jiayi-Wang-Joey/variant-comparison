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
    "isoLASER"     = "#FDBF6F"
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
        strip.text = element_markdown(size = 9),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    )
)

# ---- shared inputs ----
redi_files <- args[[1]][!grepl("/spectrum/", args[[1]])]
redi <- rbindlist(lapply(redi_files, fread))
redi <- redi[Type == "SNP" & tool != "GATK"]
redi[, coverage := as.integer(coverage)]
redi <- add_tech(redi)

spec_files <- args[[1]][grepl("/spectrum/", args[[1]])]
spec <- rbindlist(lapply(spec_files, fread))
spec <- spec[tool != "GATK" & !grepl(",", substitution)]
spec <- add_tech(spec)

# NOTE: not produced by any Snakemake rule (untracked, manually cached RDS) -
# regenerate by hand if results/happy/ changes.
happy <- as.data.table(readRDS("data/results/happy.csv"))
happy <- happy[Filter == "PASS" & Type == "SNP" & tool != "GATK"]
happy[, coverage := as.integer(as.character(coverage))]
happy <- add_tech(happy)

# --- a: total A>G/T>C FP, stacked (REDIportal-confirmed vs other), coverage=5, by chemistry ---
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
              inherit.aes = FALSE, vjust = -0.3, size = 2) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    scale_fill_manual(values = cat_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    common_theme +
    labs(x = "Variant Caller", y = "Number of A>G/T>C FP SNVs", fill = "Category")

# --- b: improved F1 (masked - original), vs coverage, by chemistry ---
db <- merge(happy, redi[, .(tool, bamtype, aligner, sample, coverage, FP_in_REDIportal_strict)],
            by = c("tool", "bamtype", "aligner", "sample", "coverage"), all.x = TRUE)
db[is.na(FP_in_REDIportal_strict), FP_in_REDIportal_strict := 0]

db_agg <- db[, .(TP = sum(TRUTH.TP), FP = sum(QUERY.FP), FN = sum(TRUTH.FN),
                  FP_redi_strict = sum(FP_in_REDIportal_strict)),
             by = .(tool, cell_line, tech, coverage)]
db_agg[, masked_FP := FP - FP_redi_strict]
db_agg[, orig_precision := TP / (TP + FP)]
db_agg[, masked_precision := TP / (TP + masked_FP)]
db_agg[, improved_precision := masked_precision - orig_precision]
db_agg[, coverage := factor(coverage, levels = sort(unique(coverage)))]

pb <- ggplot(db_agg, aes(coverage, improved_precision, color = tool, group = tool)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    scale_color_manual(values = tool_cols) +
    common_theme +
    theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5)) +
    labs(x = "Coverage Cutoff (DP >= n)", y = "Improved Precision", color = "Variant Caller")

gg <- pa + pb + plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width = 27, height = 20, units = "cm")
