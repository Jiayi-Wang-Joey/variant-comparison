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
dt <- dt[Filter=="PASS" & coverage==5 & QQ !="*"]
dt <- dt[Subtype == "*"]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[, QQ := as.numeric(QQ)]
dt <- dt[order(QQ), by=.(method, Type, sample)]


gg <- ggplot(dt, 
             aes(x = METRIC.Recall, 
                 y = METRIC.Precision,
                 group = method)) +  
    geom_line(aes(linetype = bamtype), col="grey") + 
    facet_grid2(Type ~ sample, scales = "free", axes="all", independent="all") +
    scale_color_brewer(palette = "Paired") + 
    labs(
        title = "Precision-Recall by Quality Score",
        x = "Recall",
        y = "Precision",
        color = "Variant Caller",
        shape = "Aligner",
        linetype = "Bam type"
    ) +
    theme_minimal() +
    theme(
        # legend.position = "bottom",      
        # legend.direction = "horizontal",
        panel.border = element_rect(fill=NA),
        axis.text.x = element_text(size=8, color="grey50"), 
        axis.text.y = element_text(size=8, color="grey50")   
    )

nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.5*nSample, height=24, units="cm")
