suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(ggtext)
})
setwd("/data/jiayiwang/variant-comparison")

# y-positions for stacked-bar count labels, nudged apart to avoid overlap
add_label_pos <- function(dt, by, min_gap_frac = 0.18) {
    dt <- copy(dt)
    setorderv(dt, c(by, "Category"))
    dt[, ymax := rev(cumsum(rev(Count))), by = by]
    dt[, ymin := ymax - Count]
    dt[, ymid := (ymin + ymax) / 2]
    dt[, total := sum(Count), by = by]
    dt[, label_y := {
        o <- order(ymid)
        y <- ymid[o]
        mg <- total[1] * min_gap_frac
        n <- length(y)
        if (n > 1) {
            for (i in 2:n) if (y[i] - y[i - 1] < mg) y[i] <- y[i - 1] + mg
        }
        out <- numeric(n)
        out[o] <- y
        out
    }, by = by]
    dt
}

tools <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn", "NanoTS")
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
dt <- dt[!(grepl("longcallR|NanoTS", tool) & Type == "INDEL")]
# isoLASER indels: PacBio only
dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & sample %in% c("cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"))]
dt[sample == "Baylor-IsoSeq", sample := "IsoSeq"]  # display name
dt[, sample := factor(sample, levels = sub("Baylor-IsoSeq", "IsoSeq", samples))]
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
    counts_lab <- add_label_pos(td[Type == type_lab], by = c("tool", "aligner", "sample"))
    stack_top <- counts_lab[, .(top_y = max(label_y), total = total[1]),
                             by = .(tool, aligner, sample)]
    f1 <- merge(
        totals[Type == type_lab, .(tool, aligner, sample, F1)],
        stack_top, by = c("tool", "aligner", "sample")
    )
    f1_offset <- if (type_lab == "SNP") 8000 else 2000
    f1[, label_y := top_y + f1_offset]
    ggplot() +
        geom_col(
            data = td[Type == type_lab],
            aes(aligner, Count, fill = Category)
        ) +
        geom_text(
            data = counts_lab,
            aes(aligner, label_y, label = Count),
            size = 3.2
        ) +
        geom_text(
            data = f1,
            aes(aligner, label_y, label = sprintf("%.2f", F1)),
            size = 3.8,
            fontface = "bold"
        ) +
        facet_grid2(tool ~ sample, scales = "free") +
        scale_fill_brewer(palette = "Paired") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
        theme_classic(base_size = 14) +
        theme(
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 11),
            axis.title = element_text(size = 13),
            plot.title = element_text(size = 15, face = "bold"),
            strip.text = element_markdown(size = 12),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 11)
        ) +
        labs(y = paste0("Count (", disp_lab, ", PASS, coverage=5)"), x = "Aligner", fill = "Category",
             title = paste0(disp_lab, ": minimap2 vs minimap2-nojunc (HG004)"))
}

p_snp <- .mkplot("SNP")
p_indel <- .mkplot("INDEL")

pdf(args[[2]], width = 34 / 2.54, height = 27 / 2.54)
print(p_snp)
print(p_indel)
dev.off()
cat("done\n")
