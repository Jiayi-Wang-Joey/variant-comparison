suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
.f <- \(filter, title) {
    dt <- dt[Filter==filter & 
                 coverage %in% c("gt5", "gt10", "gt30", "gt50", "gt100")]
    dt <- dt[!is.na(METRIC.F1_Score)]
    dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
    dt[,coverage:=substr(coverage, 3, nchar(coverage))]
    dt[,coverage:=factor(coverage, levels = c("5", "10", "30", "50", "100"))]
    dt[,method:=paste(bamtype, aligner, tool, sep=">")]
    
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
            x = "Coverage Cutoff (DP >= n)",
            y = "F1 Score",
            color = "Method",
        ) +
        scale_color_brewer(palette = "Set2") +
        theme(
            panel.border=element_rect(fill=NA)) +
        ggtitle(title)
}

p1 <- .f("PASS", "F1 score - PASS only")
p2 <- .f("ALL", "F1 score - ALL")
gg <- (p1 / p2) + plot_layout(guides = "collect")
nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.5*nSample, height=24, units="cm")