# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.roc.csv.gz$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & coverage=="gt5" & QQ !="*"]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[, QQ := as.numeric(QQ)]
dt <- dt[order(QQ), by=.(method, Type, sample)]


gg <- ggplot(dt,
             aes(x = METRIC.Recall, 
                 y = METRIC.Precision,
                 color = method)) +
    geom_line(alpha = 0.5) +
    facet_grid2(Type ~ sample, scales="free", axes="all") +
    scale_color_brewer(palette = "Paired") +
    labs(
        title = "Precision-Recall by Sample and Variant Type (PASS-only)",
        x = "Recall",
        y = "Precision",
        color = "Method"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        panel.border=element_rect(fill=NA)) +
    guides(
        color = guide_legend(nrow = 4, override.aes = list(size = 3)),
        size = guide_legend(nrow = 3)
    )  

ggsave(args[[2]], gg, width=25, height=20, units="cm")
