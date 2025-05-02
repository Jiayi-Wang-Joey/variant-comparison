# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/stratified/",
# pattern = "\\.extended.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/stratified-GCcontent.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(dplyr)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS" & 
             coverage %in% c("gt5", "gt10", "gt30", "gt50", "gt100")]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt <- dt[coverage=="gt5"]
dt <- dt[grepl("^gc[0-9]", Subset)]
dt[,method:=paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(dt$Subset,
                    "gc15_slop50"     = "<15",
                    "gc15to20_slop50" = "15-20",
                    "gc20to25_slop50" = "20-25",
                    "gc25to30_slop50" = "25-30",
                    "gc30to55_slop50" = "30-55",
                    "gc55to60_slop50" = "55-60",
                    "gc60to65_slop50" = "60-65",
                    "gc65to70_slop50" = "65-70",
                    "gc70to75_slop50" = "70-75",
                    "gc75to80_slop50" = "75-80",
                    "gc80to85_slop50" = "80-85",
                    "gc85_slop50"     = ">85"
)

lvls <- c("<15", "15-20", "20-25", "25-30", "30-55", "55-60",
    "60-65", "65-70", "70-75", "75-80", "80-85", ">85")
dt$Subset <- factor(dt$Subset, levels = lvls)

shown <- lvls[c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, 
                FALSE, TRUE, FALSE, TRUE, TRUE)]


gg <- ggplot(dt, aes(Subset, METRIC.F1_Score, 
                     col=tool,
                     linetype = aligner,
                     shape = bamtype,
                     group = method)) +
    geom_point(alpha=0.6) +
    geom_line(alpha=0.6) + 
    facet_grid2(Type ~ sample, scales = "free") +
    theme_minimal() +
    scale_x_discrete(breaks = shown) +
    labs(
        x = "GC content %",
        y = "F1 Score",
        color = "Method",
    ) +
    scale_color_brewer(palette = "Set2") +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border=element_rect(fill=NA))
nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=5.7*nSample, height=14, units="cm")