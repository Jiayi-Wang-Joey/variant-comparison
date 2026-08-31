suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})
setwd("/data/jiayiwang/variant-comparison")

aligners <- c("minimap2", "minimap2-nojunc", "pbmm2")
samples <- c("MasSeq", "Baylor-IsoSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004")

rows <- list()
for (aln in aligners) {
    for (s in samples) {
        f <- file.path("results/align/raw/supplementary_stats", paste0(aln, "_HG004-", s, ".tsv"))
        if (file.exists(f)) {
            dt <- fread(f)
            val <- dt[metric == "supplementary_read_rate", value]
            if (length(val) == 1) {
                rows[[length(rows) + 1]] <- data.table(aligner = aln, sample = s, supplementary_read_rate = val)
            }
        }
    }
}
dt <- rbindlist(rows)
dt[, sample := factor(sample, levels = samples)]
dt[, aligner := factor(aligner, levels = aligners)]

cols <- c("minimap2" = "#d95f02", "minimap2-nojunc" = "#7570b3", "pbmm2" = "#1b9e77")

gg <- ggplot(dt, aes(x = sample, y = supplementary_read_rate, fill = aligner)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = round(supplementary_read_rate, 2)),
              position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
    scale_fill_manual(values = cols) +
    labs(x = "HG004 Dataset", y = "Supplementary Read Rate (%)",
         fill = "Aligner", title = "Supplementary Read Rate by Aligner (HG004, raw)") +
    theme_classic() +
    theme(
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
    )

ggsave("plts/supp_read_rate_HG004_by_aligner.pdf", gg, width = 20, height = 12, units = "cm")
fwrite(dt, "data/results/supp_read_rate_HG004_by_aligner.csv")
cat("done\n")
