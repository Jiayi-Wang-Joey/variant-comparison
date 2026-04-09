suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(RColorBrewer)
    library(tidytext)
})

dt <- fread("results/stratified/GATK_origin_minimap2_HG002-dRNA002_1.extended.csv")
dt <- unique(dt[,.(Subset, Subtype, Subset.Size)])
dt <- dt[Subtype=="*"]
homo <- dt[grepl("homopolymer", Subset)]
homo$Subset <- recode(homo$Subset,
                      # "SimpleRepeat_homopolymer_4to6_AT_slop5"   = "4–6 AT",
                      # "SimpleRepeat_homopolymer_4to6_GC_slop5"   = "4–6 GC",
                      "SimpleRepeat_homopolymer_4to6_slop5"      = "4–6",
                      
                      # "SimpleRepeat_homopolymer_7to11_AT_slop5"  = "7–11 AT",
                      # "SimpleRepeat_homopolymer_7to11_GC_slop5"  = "7–11 GC",
                      "SimpleRepeat_homopolymer_7to11_slop5"     = "7–11",
                      
                      # "SimpleRepeat_homopolymer_ge12_AT_slop5"   = "≥12 AT",
                      # "SimpleRepeat_homopolymer_ge12_GC_slop5"   = "≥12 GC",
                      "SimpleRepeat_homopolymer_ge12_slop5"      = "≥12",
                      
                      # "SimpleRepeat_homopolymer_ge21_AT_slop5"   = "≥21 AT",
                      # "SimpleRepeat_homopolymer_ge21_GC_slop5"   = "≥21 GC",
                      "SimpleRepeat_homopolymer_ge21_slop5"      = "≥21"
)
homo$Subset <- factor(homo$Subset,
                levels = c("4–6", "7–11", "≥12", "≥21"
                           #"4–6 AT", "7–11 AT", "≥12 AT",  "≥21 AT", 
                           #"4–6 GC", "7–11 GC", "≥12 GC", "≥21 GC"
                          ))
homo <- na.omit(homo)

keep_subsets <- c(
    #"AllTandemRepeatsandHomopolymers_slop5",
    # "AllTandemRepeats_lt51bp_slop5",
    # "AllTandemRepeats_51to200bp_slop5",
    # "AllTandemRepeats_201to10000bp_slop5",
    # "AllTandemRepeats_gt10000bp_slop5"
    "SimpleRepeat_diTR_10to49_slop5",
    "SimpleRepeat_diTR_50to149_slop5",
    "SimpleRepeat_triTR_14to49_slop5",	
    "SimpleRepeat_triTR_50to149_slop5",	
    "SimpleRepeat_quadTR_19to49_slop5", 
    "SimpleRepeat_quadTR_50to149_slop5"
)

tdm <- dt[Subset %in% keep_subsets]

tdm$Subset <- recode(tdm$Subset,
                     "SimpleRepeat_diTR_10to49_slop5"="diTR_10-49",
                     "SimpleRepeat_diTR_50to149_slop5"="diTR_50-149",
                     "SimpleRepeat_triTR_14to49_slop5"="triTR_14-49",	
                     "SimpleRepeat_triTR_50to149_slop5"="triTR_50-149",	
                     "SimpleRepeat_quadTR_19to49_slop5"="quadTR_19-49", 
                     "SimpleRepeat_quadTR_50to149_slop5"="quadTR_50-149",
                    # "AllTandemRepeats_lt51bp_slop5" = "<51",
                    # "AllTandemRepeats_51to200bp_slop5" = "51-200",
                    # "AllTandemRepeats_201to10000bp_slop5" = "201-10000",
                    # "AllTandemRepeats_gt10000bp_slop5" = ">10000"
)

lvls <- c(
    "diTR_10-49", "diTR_50-149",
    "triTR_14-49", "triTR_50-149",
    "quadTR_19-49", "quadTR_50-149"
    #"<51", "51-200", "201-10000", ">10000"
)

tdm$Subset <- factor(tdm$Subset, levels = lvls)

tdm$Subset <- factor(tdm$Subset, levels = lvls)
gc <- dt[grepl("^gc[0-9]", Subset)]
gc$Subset <- recode(gc$Subset,
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
gc$Subset <- factor(gc$Subset, levels = lvls)

homo$type <- "Homopolymer" 
tdm$type <- "Tandeam Repeat" 
gc$type <- "GC Content"

td <- do.call(rbind, list(homo, tdm, gc))
td$type <- factor(td$type, levels = c("Homopolymer", 
                                      "Tandeam Repeat", 
                                      "GC Content"))
cols <- colorRampPalette(brewer.pal(12, "Paired"))(3)
gg <- ggplot(td, aes(Subset, Subset.Size, fill = type)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~type, scales = "free", ncol=3) +
    theme_classic() +
    scale_color_manual(values = cols) +
    labs(
        x = "Homopolymer / Tandem Repeats length or GC content (%)",
        y = "Genomic size (bp)"
    ) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
    
    
ggsave("plts/stratified_region.pdf", gg, height = 8, width=18, units = "cm")


