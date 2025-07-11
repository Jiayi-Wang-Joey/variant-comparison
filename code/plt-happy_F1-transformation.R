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
dt <- dt[coverage==5,]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,Method:=paste(aligner,tool,sep=".")]
dt[,bamtype:=factor(bamtype, levels = c("origin", "transformed"))]

gg <- ggplot(dt, aes(x = bamtype, y = METRIC.F1_Score, color = tool, 
                          group = Method)) +
    geom_line() + geom_point(aes(shape=aligner)) +
    facet_grid2(Type ~ sample, scales = "free_y", axes = "y", 
                independent = "y") +  
    scale_color_brewer(palette = "Set1") + 
    theme_minimal() +
    labs(
        x = "Bam Type",
        y = "F1 Score",
        color = "Caller",
        shape = "Aligner"
    ) +
    theme(
        legend.position = "bottom",      
        legend.direction = "horizontal",
        panel.border=element_rect(fill=NA),
        strip.text = element_text(size = 10)   
    )

nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.2*nSample, height=12, units="cm")