# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]


dt[,method:=paste(bamtype, aligner, tool, sep=">")]


gg <- ggplot(dt, aes(reorder(method,METRIC.F1_Score), 
                     METRIC.F1_Score, fill = method)) +
    geom_boxplot(alpha=0.6) +
    facet_grid2(Type ~ sample, scales = "free") +
    theme_minimal() +
    labs(
        title = "Comparison of Methods F1 Scores",
        x = "Method",
        y = "F1 Score",
        fill = "Variant caller",
    ) +
    scale_fill_brewer(palette = "Paired") +
    #guides(fill = guide_legend(nrow = 4, override.aes = list(size = 3))) +
    theme(
          axis.text.x = element_blank(),
          panel.border=element_rect(fill=NA))

ggsave(args[[2]], gg, width=30, height=16, units="cm")