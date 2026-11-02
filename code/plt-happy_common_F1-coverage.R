# standard vs. common evaluation region, HG004, minimap2
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(ggtext)
    library(scales)
    library(patchwork)
})

fs <- args[[1]]
res <- lapply(fs, \(f) {
    x <- fread(f, header = TRUE)
    if (nrow(x)) x[, region := ifelse(grepl("happy_common", f), "Common", "Standard")]
    x
})
dt <- rbindlist(res[vapply(res, nrow, integer(1)) > 0])
dt <- dt[Filter == "PASS" & !is.na(METRIC.F1_Score)]
dt <- dt[aligner == "minimap2" & bamtype == "origin" & grepl("^HG004", sample)]
dt <- dt[!(tool %in% c("longcallR", "longcallR-nn", "NanoTS") & Type == "INDEL")]
dt <- dt[!(Type == "INDEL" & tool == "isoLASER" & grepl("dRNA|cDNA", sample))]
dt[, coverage := factor(coverage, levels = sort(unique(coverage)))]
dt[, Type := factor(Type, c("SNP", "INDEL"))]
dt[, region := factor(region, c("Standard", "Common"))]
dt[, platform := sub("^[^-]+-", "", sample)]
dt[platform == "Baylor-IsoSeq", platform := "IsoSeq"]
dt[, platform := factor(platform, c("IsoSeq", "MasSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"))]
fwrite(dcast(dt, tool + sample + Type + coverage ~ region,
             value.var = c("METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall", "TRUTH.TOTAL")),
       "data/results/happy_common_metrics.csv")

cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F",
    "NanoTS"       = "#E05A8A"
)
metric_levels <- c("F1", "Precision", "Recall", "Truth variants")

long <- melt(dt, id.vars = c("tool", "platform", "Type", "coverage", "region"),
             measure.vars = c("METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall"),
             variable.name = "metric", value.name = "value")
long[, metric := factor(sub("METRIC\\.", "", sub("_Score", "", metric)), metric_levels)]
n_truth <- dt[, .(value = max(TRUTH.TOTAL)), by = .(platform, Type, coverage, region)]
n_truth[, metric := factor("Truth variants", metric_levels)]

.page <- \(type) {
    ggplot(long[Type == type], aes(coverage, value, linetype = region,
                                   group = interaction(tool, region))) +
        geom_line(aes(colour = tool), linewidth = 0.5) +
        geom_point(aes(colour = tool), size = 1) +
        geom_line(data = n_truth[Type == type], aes(coverage, value, linetype = region, group = region),
                  inherit.aes = FALSE, colour = "grey30", linewidth = 0.5) +
        geom_point(data = n_truth[Type == type], aes(coverage, value),
                   inherit.aes = FALSE, colour = "grey30", size = 1) +
        facet_grid2(metric ~ platform, scales = "free_y",
                    labeller = labeller(platform = \(x) ifelse(x %in% c("IsoSeq", "MasSeq"),
                        paste0("<span style='color:#54278f;'>", x, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", x, "</span>")))) +
        facetted_pos_scales(y = list(
            metric %in% c("F1", "Precision", "Recall") ~ scale_y_continuous(limits = c(0, 1)),
            metric == "Truth variants" ~ scale_y_continuous(labels = label_comma())
        )) +
        scale_colour_manual(values = cols) +
        scale_linetype_manual(values = c(Standard = "solid", Common = "dashed")) +
        theme_bw(base_size = 11) +
        theme(strip.text = element_markdown(),
              strip.background = element_rect(fill = "white"),
              panel.grid.minor = element_blank()) +
        labs(x = "Coverage Cutoff (DP >= n)", y = NULL, colour = "Variant Caller", linetype = "Region",
             title = ifelse(type == "SNP", "SNV", "INDEL"))
}

p <- .page("SNP") / (.page("INDEL") + theme(legend.position = "none")) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 16))

pdf(args[[2]], width = 32 / 2.54, height = 36 / 2.54)
print(p)
dev.off()
