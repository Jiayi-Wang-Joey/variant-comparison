suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(RColorBrewer)
    library(ggtext)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[
    coverage == 5 &
        Filter == "PASS" &
        !is.na(METRIC.F1_Score) &
        !(grepl("longcallR", tool) & Type == "INDEL") &
        grepl("MasSeq|IsoSeq", sample)
]
dt[,FDR:=QUERY.FP/(TRUTH.TP+QUERY.FP)]
dt[,FNR:=TRUTH.FN/(TRUTH.TP+TRUTH.FN)]
m <- melt(
    dt,
    id.vars = setdiff(names(dt), c("FDR","FNR")),
    measure.vars = c("FDR","FNR"),
    variable.name = "metric",
    value.name = "val"
)

td <- data.table::dcast(m, sample + tool + Type + metric ~ aligner, 
                        value.var = "val")
td[,Type:=factor(Type, levels=c("SNP", "INDEL"))]
saveRDS(td, "data/results/aligner.rds")
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER" = "#FDBF6F"
)
td[, sample:=gsub("-Baylor", "", sample)]
td[grepl("HG002a", sample), sample:=sub("HG002a", "HG002", sample)]
gg <- ggplot(td, aes(minimap2, pbmm2, color=tool, shape=metric)) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(alpha=0.8, size=3) +
    facet_grid(Type~sample) +
    labs(
        x = "minimap2",
        y = "pbmm2",
        color = "Variant Caller",
        shape = "Metric"
    ) +
    theme_classic() +
    scale_color_manual(values = cols) +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.6
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.text  = element_text(size = 9)) + coord_equal()

ggsave(args[[2]], gg, width=30, height=12, units="cm")