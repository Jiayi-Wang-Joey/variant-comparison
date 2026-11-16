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
    dt <- dt[cell_line %in% c("HG002", "HG004", "HG005")]
    dt[, cell_line := factor(cell_line, levels = c("HG002", "HG004", "HG005"))]
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

# change in precision and recall after removing calls at REDIportal sites
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
d[, F1_masked := 2 * TP_masked / (2 * TP_masked + FP_masked + FN_masked)]
fwrite(d[, tech := gsub("<[^>]+>", "", tech)], sub("\\.pdf$", ".csv", args[[2]]))
d <- add_tech(d[, sample := paste0(cell_line, "-", tech)])
d[, coverage := factor(coverage, levels = sort(unique(coverage)))]

.p <- \(y, lab, dd = d, zero = TRUE) ggplot(dd, aes(coverage, .data[[y]], color = tool, group = tool)) +
    (if (zero) geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    scale_color_manual(values = tool_cols) +
    common_theme +
    theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) +
    labs(x = "Coverage Cutoff (DP >= n)", y = lab, color = "Variant Caller")

pc <- .p("F1_masked", "F1 (REDIportal-filtered)", zero = FALSE)
gg <- .p("dPrecision", expression(atop(Delta * " Precision", "(REDIportal-filtered)")), d[cell_line == "HG002"]) /
      .p("dRecall", expression(atop(Delta * " Recall", "(REDIportal-filtered)")), d[cell_line == "HG002"]) /
      pc +
    plot_layout(guides = "collect", heights = c(1, 1, 3)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 20, face = "bold"), legend.position = "bottom")

ggsave(args[[2]], gg, width = 28, height = 32, units = "cm")
