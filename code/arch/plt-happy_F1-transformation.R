suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]

gg <- ggplot(dt, aes(reorder(tool,METRIC.F1_Score), 
                     METRIC.F1_Score, 
                     fill = reorder(bamtype, METRIC.F1_Score))) +
    geom_boxplot(position=position_dodge(0.8), alpha=0.6) +
    facet_grid2(Type ~ sample, scales = "free_y", axes = "y") +
    theme_bw() +
    labs(
        title = "Comparison of Bam Transformation: F1 Scores",
        x = "Variant caller",
        y = "F1 Score",
        fill = "Bam Type"
    ) +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(args[[2]], gg, width=26, height=16, units="cm")