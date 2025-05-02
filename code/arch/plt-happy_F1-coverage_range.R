suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & coverage %in% c("1-5", "5-20", "20-50", 
                                          "50-100", "gt100")]
dt <- dt[!is.na(METRIC.F1_Score)]
dt[,coverage:=factor(coverage, levels = c("1-5", "5-20", "20-50", 
                                          "50-100", "gt100"))]

dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]

gg <- ggplot(dt, aes(coverage, 
                     METRIC.F1_Score, 
                     color = tool,
                     linetype = aligner,
                     shape = bamtype,
                     group = method)) +
    geom_point(alpha=0.6) +
    geom_line(alpha=0.6) + 
    facet_grid2(Type ~ sample, scales = "free") +
    theme_minimal() +
    labs(
        x = "coverage",
        y = "F1 Score",
        color = "Method",
    ) +
    scale_color_brewer(palette = "Set2") +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border=element_rect(fill=NA))

ggsave(args[[2]], gg, width=30, height=14, units="cm")