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
dt <- dt[Filter=="PASS" & 
             coverage %in% c("gt5", "gt10", "gt30", "gt50", "gt100")]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,coverage:=substr(coverage, 3, nchar(coverage))]
dt[,coverage:=factor(coverage, levels = c("5", "10", "30", "50", "100"))]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[,aligner_tool:=paste(aligner, tool, sep = ".")]

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
        color = "Coverage cutoff"
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
ggsave(args[[2]], gg, width=5*nSample, height=12, units="cm")
