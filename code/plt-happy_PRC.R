# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,coverage:=factor(coverage)]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[,aligner_tool:=paste(aligner, tool, sep = ".")]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
sample_levels <- dt[, .(platform = unique(platform)), by = sample][
    order(platform)
]$sample
dt[, sample := factor(sample, levels = sample_levels)]

coverage_colors <- c(
    "5" = "#4575b4",
    "10" = "#91bfdb",
    "30" = "#fee090",
    "50" = "#fc8d59",
    "100" = "#d73027"
)

gg <- ggplot(dt, 
       aes(x = METRIC.Recall, 
           y = METRIC.Precision,
           group = method)) +  
    geom_line(aes(linetype = bamtype), col="grey") + 
    geom_point(aes(shape = aligner_tool, col = coverage)) +  
    facet_grid2(Type ~ sample, scales = "free", axes="all", independent="all") +
    scale_color_brewer(palette = "Paired") + 
    scale_color_manual(values = coverage_colors) +  
    labs(
        title = "Precision-Recall by Sample and Variant Type (PASS-only)",
        x = "Recall",
        y = "Precision",
        color = "Coverage cutoff",
        shape = "Aligner.Caller",
        linetype = "Bam type"
    ) +
    theme_minimal() +
    scale_x_continuous(
        labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(
        labels = label_number(accuracy = 0.01)) +
    theme(
        legend.position = "bottom",      
        legend.direction = "horizontal",
        panel.border = element_rect(fill=NA),
        axis.text.x = element_text(size=8, color="grey50"), 
        axis.text.y = element_text(size=8, color="grey50")   
    )

                
nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.5*nSample, height=14, units="cm")
