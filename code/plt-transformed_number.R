suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(ggtext)
    library(patchwork)
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
        !is.na(METRIC.F1_Score) &
        aligner == "minimap2" &
        !(grepl("longcallR|NanoTS", tool) & Type == "INDEL")
]
dt[, platform := sapply(strsplit(sample, "-"), tail, 1)]
dt[, platform := factor(platform, levels = c("IsoSeq", "MasSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"))]
levels(dt$platform) <- ifelse(levels(dt$platform) %in% c("IsoSeq", "MasSeq"),
    paste0("<span style='color:#54278f;'>", levels(dt$platform), "</span>"),
    paste0("<span style='color:#d95f0e;'>", levels(dt$platform), "</span>"))
dt[, bamtype := factor(bamtype, levels = c("origin", "transformed"))]
td <- melt(
    dt,
    id.vars = c("Type", "tool", "platform", "bamtype"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP", "METRIC.F1_Score"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]

plot_counts <- function(d, title) {
    bars <- d[Category %in% c("TP", "FN", "FP")]
    lab <- add_label_pos(bars, by = c("tool", "platform", "bamtype"), scale_by = "tool")
    f1 <- merge(
        d[Category == "F1_Score", .(tool, platform, bamtype, F1 = Count)],
        bars[, .(Total = sum(Count)), by = .(tool, platform, bamtype)],
        by = c("tool", "platform", "bamtype")
    )
    f1[, offset := 0.05 * max(Total), by = tool]
    ggplot() +
        geom_col(data = bars, aes(bamtype, Count, fill = Category)) +
        geom_text(data = lab, aes(bamtype, label_y, label = Count), size = 2.4) +
        geom_text(
            data = f1,
            aes(bamtype, Total + offset, label = sprintf("%.2f", F1)),
            size = 2.6,
            fontface = "bold"
        ) +
        facet_grid2(tool ~ platform, scales = "free") +
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
        labs(y = "Total", x = "BAM type") +
        ggtitle(title)
}

p1 <- plot_counts(td[Type == "SNP"], "SNV")
p2 <- plot_counts(td[Type == "INDEL"], "INDEL")

gg <- p1 + p2 + plot_layout(ncol = 1, heights = c(4, 2), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width = 24, height = 34, units = "cm")
