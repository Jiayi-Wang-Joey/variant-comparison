suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})
setwd("/data/jiayiwang/variant-comparison")

tools <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn")
samples <- c("MasSeq", "Baylor-IsoSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004")
aligners <- c("minimap2", "minimap2-nojunc", "pbmm2")

rows <- list()
for (t in tools) {
    for (a in aligners) {
        for (s in samples) {
            f <- file.path("results/happy", paste0(t, "_origin_", a, "_HG004-", s, "_5.summary.csv"))
            if (file.exists(f)) {
                dt <- fread(f)
                sub <- dt[Type == "SNP" & Filter == "PASS"]
                if (nrow(sub) == 1) {
                    rows[[length(rows) + 1]] <- data.table(
                        tool = t, aligner = a, sample = s,
                        precision = sub$METRIC.Precision,
                        recall = sub$METRIC.Recall,
                        F1 = sub$METRIC.F1_Score
                    )
                }
            }
        }
    }
}
dt <- rbindlist(rows)
dt[, sample := factor(sample, levels = samples)]
dt[, aligner := factor(aligner, levels = aligners)]
dt[, tool := factor(tool, levels = tools)]

fwrite(dt, "data/results/variant_perf_HG004_by_aligner.csv")

cols <- c("minimap2" = "#d95f02", "minimap2-nojunc" = "#7570b3", "pbmm2" = "#1b9e77")

gg <- ggplot(dt, aes(x = sample, y = F1, color = aligner, group = aligner)) +
    geom_point(size = 2, alpha = 0.85) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    facet_wrap(~tool, ncol = 5) +
    scale_color_manual(values = cols) +
    labs(x = "HG004 Dataset", y = "F1 Score (SNP, PASS, coverage=5)",
         color = "Aligner", title = "Variant Calling Performance by Aligner (HG004)") +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        legend.position = "top"
    )

ggsave("plts/variant_perf_HG004_by_aligner.pdf", gg, width = 32, height = 10, units = "cm")
cat("done\n")
