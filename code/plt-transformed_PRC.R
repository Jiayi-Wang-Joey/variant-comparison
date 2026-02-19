suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(RColorBrewer)
    library(ggtext)
    library(patchwork)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[coverage==5]
nk <- 5
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk),
                 unique(dt$tool))
dt <- dt[!(grepl("longcallR", tool) & Type=="INDEL"),]

dt[,bamtype:=factor(bamtype, levels = c("origin", "transformed"))]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]

setnames(dt, "METRIC.Recall", "Recall")
setnames(dt, "METRIC.Precision", "Precision")
dt$METRIC.F1_score <- NULL

dt <- melt(
    dt,
    measure.vars = c("Recall", "Precision"),
    variable.name = "Metric",
    value.name = "Value"
)
dt[,Method:=paste(aligner,tool,Metric,sep=".")]


gg <- ggplot(dt, aes(x = bamtype,
                     y = Value,
                     color = tool,
                     shape = Metric,
                     group = Method)) +
    geom_line() +
    geom_point() +
    facet_grid2(Type ~ platform, scales = "free_y", axes = "y",
                independent = "y") +
    scale_color_manual(values = cols) +
    theme_classic() +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_markdown(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing.y = unit(0, "lines")
    ) +
    labs(x="BAM Type", color="Variant Caller")



ggsave(args[[2]], gg, width=32, height=12, units="cm")
