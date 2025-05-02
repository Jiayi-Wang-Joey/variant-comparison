suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & coverage=="gt5"]
dt <- dt[,Method1:=paste(bamtype,tool,sep=">")]
dt <- dt[,Method:=paste(bamtype,aligner,tool,sep=">")]


dt_long <- melt(dt, id.vars = c("sample", "Type", "aligner", 
                                "Method", "Method1"), 
                measure.vars = c("TRUTH.FN", "QUERY.FP", "TRUTH.TP"),
                variable.name = "Metric", 
                value.name = "Count")

gg <- ggplot(dt_long, aes(x = aligner, y = Count, color = Metric, 
                          group = interaction(Method1, Metric))) +
    geom_line() + geom_point() +
    facet_grid2(Type ~ sample, scales = "free_y", axes = "y") +  
    scale_color_brewer(palette = "Set2") + 
    theme_minimal() +
    labs(
        x = "Aligner",
        y = "Count",
        color = "Metric"
    ) +
    theme(
        panel.border=element_rect(fill=NA),
        strip.text = element_text(size = 10, face = "bold"),  
        axis.text.x = element_text(angle = 45, hjust = 1)    
    )

ggsave(args[[2]], gg, width=25, height=20, units="cm")