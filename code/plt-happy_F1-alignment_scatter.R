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
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,coverage:=factor(coverage, levels = c("5", "10", "30", "50", "100"))]
dt <- dt[grepl("IsoSeq|MasSeq", sample, ignore.case = TRUE)]
dt <- data.table::dcast(dt, sample + coverage + tool + bamtype + Type ~ aligner, 
            value.var = "METRIC.F1_Score")
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
sample_levels <- dt[, .(platform = unique(platform)), by = sample][
    order(platform)
]$sample
dt[, sample := factor(sample, levels = sample_levels)]

gg <- ggplot(dt, aes(minimap2, pbmm2, color=tool, 
                     shape=bamtype)) + 
    geom_point(alpha=0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    facet_grid(Type~sample,  scales = "free") +
    labs(
        title = "F1 Score: pbmm2 vs minimap2",
        x = "minimap2 F1 Score",
        y = "pbmm2 F1 Score"
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set2") +
    theme(
        panel.border=element_rect(fill=NA)) 

nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.5*nSample, height=12, units="cm")