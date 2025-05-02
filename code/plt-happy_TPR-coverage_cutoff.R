suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & 
             coverage %in% c("gt5", "gt10", "gt30", "gt50", "gt100")]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,coverage:=substr(coverage, 3, nchar(coverage))]
dt[,coverage:=factor(coverage, levels = c("5", "10", "30", "50", "100"))]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[,TPR:=TRUTH.TP / TRUTH.TOTAL]

gg <- ggplot(dt, aes(coverage, 
                     TPR, 
                     color = tool,
                     linetype = aligner,
                     shape = bamtype,
                     group = method)) +
    geom_point(alpha=0.6) +
    geom_line(alpha=0.6) + 
    facet_grid2(Type ~ sample, scales = "free") +
    theme_minimal() +
    labs(
        x = "Coverage Cutoff (DP >= n)",
        y = "TPR",
        color = "Method",
    ) +
    scale_color_brewer(palette = "Set2") +
    theme(
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border=element_rect(fill=NA))
nSample <- length(unique(dt$sample))

ggsave(args[[2]], gg, width=5.5*nSample, height=12, units="cm")