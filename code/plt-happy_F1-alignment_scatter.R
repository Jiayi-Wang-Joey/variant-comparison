suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(RColorBrewer)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[coverage==5]
# dt[,coverage:=factor(coverage, levels = c("5", "10", "30", "50", "100"))]
dt <- dt[grepl("IsoSeq|MasSeq", sample, ignore.case = TRUE)]
dt <- data.table::dcast(dt, sample + coverage + tool + bamtype + Type ~ aligner, 
            value.var = "METRIC.F1_Score")
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
sample_levels <- dt[, .(platform = unique(platform)), by = sample][
    order(platform)
]$sample
dt[, sample := factor(sample, levels = sample_levels)]
dt[, Type:=factor(Type, levels = "SNP", "INDEL")]
nk <- length(unique(dt$tool))
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk),
                 unique(dt$tool))
gg <- ggplot(dt, aes(minimap2, pbmm2, color=tool)) + 
    geom_point(alpha=0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    facet_grid2(Type~sample,  scales = "free") +
    labs(
        x = "minimap2 F1 Score",
        y = "pbmm2 F1 Score",
        color = "Variant Caller"
    ) +
    theme_classic() +
    scale_color_manual(values = cols) +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) 


ggsave(args[[2]], gg, width=30, height=12, units="cm")