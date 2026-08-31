suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggh4x)
    library(patchwork)
    library(RColorBrewer)
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
                sub <- dt[Type == "SNP" & Filter == "PASS"]
                if (nrow(sub) == 1) {
                    rows[[length(rows) + 1]] <- data.table(
                        tool = t, aligner = a, sample = s,
                        TP = sub$TRUTH.TP, FN = sub$TRUTH.FN, FP = sub$QUERY.FP
                    )
                }
            }
        }
    }
}
dt <- rbindlist(rows)
dt[, FDR := FP / (TP + FP)]
dt[, FNR := FN / (TP + FN)]
dt[, sample := factor(sample, levels = samples)]
dt[, tool := factor(tool, levels = tools)]

fwrite(dt, "data/results/happy_TP-FN-FP_HG004_minimap2-vs-nojunc.csv")

tool_cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "isoLASER"     = "#FDBF6F",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)

common_theme <- list(
    theme_classic(),
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_text(size = 9),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        aspect.ratio = 1
    )
)

# --- panel a: FDR/FNR scatter, minimap2 vs minimap2-nojunc ---
rate_long <- melt(dt, id.vars = c("tool", "aligner", "sample"),
                   measure.vars = c("FDR", "FNR"), variable.name = "metric", value.name = "val")
rate_wide <- dcast(rate_long, tool + sample + metric ~ aligner, value.var = "val")

pa <- ggplot(rate_wide, aes(minimap2, `minimap2-nojunc`, color = tool, shape = metric)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(alpha = 0.8, size = 2.5) +
    facet_wrap(~sample, nrow = 1) +
    labs(x = "minimap2", y = "minimap2-nojunc", color = "Variant Caller", shape = "Metric",
         title = "FDR / FNR: minimap2 vs minimap2-nojunc") +
    scale_color_manual(values = tool_cols) +
    common_theme

# --- panel b: TP/FN/FP counts, minimap2 vs minimap2-nojunc ---
count_long <- melt(dt, id.vars = c("tool", "aligner", "sample"),
                    measure.vars = c("TP", "FN", "FP"), variable.name = "metric", value.name = "val")
count_wide <- dcast(count_long, tool + sample + metric ~ aligner, value.var = "val")
count_wide[, metric := factor(metric, levels = c("TP", "FN", "FP"))]

pb <- ggplot(count_wide, aes(minimap2, `minimap2-nojunc`, color = tool)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(alpha = 0.8, size = 2.5) +
    facet_grid2(metric ~ sample, scales = "free") +
    labs(x = "minimap2", y = "minimap2-nojunc", color = "Variant Caller",
         title = "TP / FN / FP counts: minimap2 vs minimap2-nojunc") +
    scale_color_manual(values = tool_cols) +
    common_theme +
    theme(aspect.ratio = 1)

gg <- pa / pb + plot_layout(heights = c(1, 2.2)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave("plts/happy_FDR-FNR_minimap2-vs-nojunc_HG004.pdf", gg, width = 34, height = 26, units = "cm")
cat("done\n")
