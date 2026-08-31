suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(dplyr)
    library(patchwork)
    library(ggtext)
    library(RColorBrewer)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[metric == "callable_bases"]

gc_lvls  <- c("<15", "15-20", "20-25", "25-30", "30-55", "55-60",
              "60-65", "65-70", "70-75", "75-80", "80-85", ">85")
homo_lvls <- c("4–6", "7–11", "≥12", "≥21")
tr_lvls  <- c("<50", "51-200", "201-10000", ">10000")

dt$category <- dplyr::recode(dt$subset,
    "gc15_slop50"     = "GC content", "gc15to20_slop50" = "GC content",
    "gc20to25_slop50" = "GC content", "gc25to30_slop50" = "GC content",
    "gc30to55_slop50" = "GC content", "gc55to60_slop50" = "GC content",
    "gc60to65_slop50" = "GC content", "gc65to70_slop50" = "GC content",
    "gc70to75_slop50" = "GC content", "gc75to80_slop50" = "GC content",
    "gc80to85_slop50" = "GC content", "gc85_slop50"     = "GC content",
    "SimpleRepeat_homopolymer_4to6_slop5"  = "Homopolymer",
    "SimpleRepeat_homopolymer_7to11_slop5" = "Homopolymer",
    "SimpleRepeat_homopolymer_ge12_slop5"  = "Homopolymer",
    "SimpleRepeat_homopolymer_ge21_slop5"  = "Homopolymer",
    "AllTandemRepeats_le50bp_slop5"        = "Tandem Repeat",
    "AllTandemRepeats_51to200bp_slop5"     = "Tandem Repeat",
    "AllTandemRepeats_201to10000bp_slop5"  = "Tandem Repeat",
    "AllTandemRepeats_ge10001bp_slop5"     = "Tandem Repeat"
)
dt$category <- factor(dt$category, levels = c("GC content", "Homopolymer", "Tandem Repeat"))

dt$bin <- dplyr::recode(dt$subset,
    "gc15_slop50"     = "<15",    "gc15to20_slop50" = "15-20",
    "gc20to25_slop50" = "20-25",  "gc25to30_slop50" = "25-30",
    "gc30to55_slop50" = "30-55",  "gc55to60_slop50" = "55-60",
    "gc60to65_slop50" = "60-65",  "gc65to70_slop50" = "65-70",
    "gc70to75_slop50" = "70-75",  "gc75to80_slop50" = "75-80",
    "gc80to85_slop50" = "80-85",  "gc85_slop50"     = ">85",
    "SimpleRepeat_homopolymer_4to6_slop5"  = "4–6",
    "SimpleRepeat_homopolymer_7to11_slop5" = "7–11",
    "SimpleRepeat_homopolymer_ge12_slop5"  = "≥12",
    "SimpleRepeat_homopolymer_ge21_slop5"  = "≥21",
    "AllTandemRepeats_le50bp_slop5"        = "<50",
    "AllTandemRepeats_51to200bp_slop5"     = "51-200",
    "AllTandemRepeats_201to10000bp_slop5"  = "201-10000",
    "AllTandemRepeats_ge10001bp_slop5"     = ">10000"
)
dt$bin <- factor(dt$bin, levels = c(gc_lvls, homo_lvls, tr_lvls))

dt <- dt[grepl("^HG004", sample)]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, value_mb := value / 1e6]

aes <- list(
    geom_bar(stat = "identity"),
    facet_grid2(category ~ platform, scales = "free_x", independent = "x"),
    theme_classic(),
    labs(
        x = "Genomic context bin",
        y = "Callable bases (Mb, DP >= 5)"
    ),
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
        strip.text = element_markdown(size = 11),
        axis.text.x = element_text(angle = 45, size = 7, hjust = 1, vjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)
    )
)

gg <- ggplot(dt, aes(bin, value_mb)) + aes

ggsave(args[[2]], gg, width=32, height=22, units="cm")
