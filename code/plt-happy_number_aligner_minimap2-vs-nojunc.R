suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(ggtext)
})
setwd("/data/jiayiwang/variant-comparison")

tools <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn")
samples <- c("MasSeq", "Baylor-IsoSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004")
aligners <- c("minimap2", "minimap2-nojunc")

rows <- list()
for (t in tools) {
    for (a in aligners) {
        for (s in samples) {
            f <- file.path("results/happy", paste0(t, "_origin_", a, "_HG004-", s, "_5.summary.csv"))
            if (file.exists(f)) {
                dt <- fread(f)
                sub <- dt[Type %in% c("SNP", "INDEL") & Filter == "PASS"]
                for (i in seq_len(nrow(sub))) {
                    rows[[length(rows) + 1]] <- data.table(
                        tool = t, aligner = a, sample = s, Type = sub$Type[i],
                        TP = sub$TRUTH.TP[i], FN = sub$TRUTH.FN[i], FP = sub$QUERY.FP[i],
                        F1 = sub$METRIC.F1_Score[i]
                    )
                }
            }
        }
    }
}
dt <- rbindlist(rows)
dt <- dt[!(grepl("longcallR", tool) & Type == "INDEL")]
# isoLASER essentially doesn't call indels on ONT chemistries (TP in the single
# digits to tens against thousands of truth indels, F1 ~0.00-0.01) -- exclude
# those panels from the INDEL page rather than show near-empty, misleading bars.
dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & sample %in% c("cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"))]
dt[, sample := factor(sample, levels = samples)]
dt[, aligner := factor(aligner, levels = aligners)]
dt[, tool := factor(tool, levels = tools)]

fwrite(dt, "data/results/happy_number-F1_HG004_minimap2-vs-nojunc.csv")

td <- melt(
    dt,
    id.vars = c("tool", "aligner", "sample", "Type", "F1"),
    measure.vars = c("TP", "FN", "FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := factor(Category, levels = c("FN", "FP", "TP"))]

totals <- td[, .(Total = sum(Count)), by = .(tool, sample, aligner, Type, F1)]

.mkplot <- \(type_lab) {
    disp_lab <- if (type_lab == "SNP") "SNV" else type_lab
    ggplot() +
        geom_col(
            data = td[Type == type_lab],
            aes(aligner, Count, fill = Category)
        ) +
        geom_text(
            data = td[Type == type_lab],
            aes(aligner, Count, label = Count, group = Category),
            position = position_stack(vjust = 0.5),
            size = 1.6
        ) +
        geom_text(
            data = totals[Type == type_lab],
            aes(aligner, Total, label = sprintf("%.2f", F1)),
            vjust = -0.4,
            size = 2.2,
            fontface = "bold"
        ) +
        facet_grid2(tool ~ sample, scales = "free") +
        scale_fill_brewer(palette = "Paired") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
        theme_classic() +
        theme(
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text = element_markdown(size = 8)
        ) +
        labs(y = paste0("Count (", disp_lab, ", PASS, coverage=5)"), x = "Aligner", fill = "Category",
             title = paste0(disp_lab, ": minimap2 vs minimap2-nojunc (HG004)"))
}

p_snp <- .mkplot("SNP")
p_indel <- .mkplot("INDEL")

pdf("plts/happy_number_minimap2-vs-nojunc_HG004.pdf", width = 32 / 2.54, height = 24 / 2.54)
print(p_snp)
print(p_indel)
dev.off()
cat("done\n")
