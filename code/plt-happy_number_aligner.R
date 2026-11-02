suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
})

# y-positions for stacked-bar count labels, nudged apart to avoid overlap
add_label_pos <- function(dt, by, scale_by = by, min_gap_frac = 0.08) {
    dt <- copy(dt)
    setorderv(dt, c(by, "Category"))
    dt[, ymax := rev(cumsum(rev(Count))), by = by]
    dt[, ymin := ymax - Count]
    dt[, ymid := (ymin + ymax) / 2]
    dt[, total := sum(Count), by = by]
    dt[, scale_ref := max(total), by = scale_by]
    dt[, label_y := {
        o <- order(ymid)
        y <- ymid[o]
        mg <- scale_ref[1] * min_gap_frac
        if (length(y) > 1) {
            for (i in 2:length(y)) if (y[i] - y[i - 1] < mg) y[i] <- y[i - 1] + mg
        }
        out <- numeric(length(y))
        out[o] <- y
        out
    }, by = by]
    dt
}

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[
    coverage == 5 &
        Filter == "PASS" &
        bamtype == "origin" &
        !is.na(METRIC.F1_Score) &
        !(grepl("longcallR|NanoTS", tool) & Type == "INDEL") &
        grepl("MasSeq|IsoSeq", sample) &
        !grepl("^HG002-", sample) &
        aligner %in% c("minimap2", "pbmm2")
]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "sample"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP", "METRIC.F1_Score"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
td[, sample:=gsub("-Baylor", "", sample)]
td[grepl("HG002a", sample), sample:=sub("HG002a", "HG002", sample)]
td <- unique(td)
saveRDS(td, "data/results/aligner.rds")
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]

snp_lab <- add_label_pos(snp[Category %in% c("TP","FN","FP")],
                          by = c("tool", "sample", "aligner"), scale_by = "tool")
idl_lab <- add_label_pos(idl[Category %in% c("TP","FN","FP")],
                          by = c("tool", "sample", "aligner"), scale_by = "tool")

p1 <- ggplot() +
    geom_col(
        data = snp[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, fill = Category)
    ) +
    geom_text(
        data = snp_lab,
        aes(aligner, label_y, label = Count),
        size = 2.4
    ) +
    geom_text(
        data = merge(
            snp[Category == "F1_Score", .(tool, sample, aligner, F1 = Count)],
            snp[Category %in% c("TP","FN","FP"),
                .(Total = sum(Count)),
                by = .(tool, sample, aligner)],
            by = c("tool","sample","aligner")
        ),
        aes(aligner, Total + 5000,
            label = sprintf("%.2f", F1)),
        size = 2.6,
        fontface = "bold"
    ) +
    facet_grid2(tool ~ sample, scales = "free") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 13) +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown(size = 10)
    ) +
    labs(y = "Total", x = "Variant caller") +
    ggtitle("SNV")

p2 <- ggplot() +
    geom_col(
        data = idl[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, fill = Category)
    ) +
    geom_text(
        data = idl_lab,
        aes(aligner, label_y, label = Count),
        size = 2.4
    ) +
    geom_text(
        data = merge(
            idl[Category == "F1_Score", .(tool, sample, aligner, F1 = Count)],
            idl[Category %in% c("TP","FN","FP"),
                .(Total = sum(Count)),
                by = .(tool, sample, aligner)],
            by = c("tool","sample","aligner")
        ),
        aes(aligner, Total + 1000,
            label = sprintf("%.2f", F1)),
        size = 2.6,
        fontface = "bold"
    ) +
    facet_grid2(tool ~ sample, scales = "free") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 13) +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown(size = 10)
    ) +
    labs(y = "Total", x = "Variant caller") +
    ggtitle("INDEL")

pdf(args[[2]], width = 24 / 2.54, height = 24 / 2.54)
print(p1)
print(p2)
dev.off()


