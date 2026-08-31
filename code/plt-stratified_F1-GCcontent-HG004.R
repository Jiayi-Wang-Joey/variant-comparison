# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/stratified/",
# pattern = "\\.extended.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/stratified-GCcontent.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(patchwork)
    library(RColorBrewer)
})

cols <- c("Filter", "METRIC.F1_Score", "TRUTH.TOTAL", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR",tool) & Type == "INDEL")]
dt <- dt[coverage==5 & Subtype == "*"]
dt <- dt[grepl("^gc[0-9]", Subset)]
dt <- dt[aligner=="minimap2"]
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
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
#dt[cell_line == "HG002a", cell_line := "HG002"]
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F"
)


dt <- dt[cell_line == "HG004"]
dt <- dt[!(Type == "INDEL" & grepl("dRNA|cDNA", platform) & tool == "isoLASER")]
dt$Type <- factor(dt$Type, levels = c("SNP", "INDEL"))
lb <- dt[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(Subset, Type, platform)]
aes <- list(geom_point(alpha=0.6),
            geom_line(alpha = 0.6),
            facet_grid(Type ~ platform),
            theme_minimal(),
            labs(x = "GC Content (%)", y = "F1 Score", color = "Variant Caller"),
            scale_color_manual(values = cols),
            guides(color = guide_legend(nrow = 1)),
            theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
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
                  strip.text =  element_markdown(size=11),
                  axis.line = element_line(color = "black", linewidth = 0.3),
                  panel.spacing = unit(0, "lines"),
                  panel.spacing.x = unit(0, "lines"),
                  panel.spacing.y = unit(0, "lines"),
                  axis.text.y = element_text(size = 7),
                  axis.title.x = element_text(size = 11),
                  axis.title.y = element_text(size = 11),
                  legend.title = element_text(size = 11),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 45, size = 7,
                                             hjust = 1, vjust = 1)))

gg <- ggplot(dt, aes(Subset, METRIC.F1_Score,
                      col=tool,
                      group = method)) + aes +
    geom_text(data = lb,
        aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
        inherit.aes = FALSE,
        angle = 30, hjust = 0.5, vjust = 1.3,
        size = 1.3, color = "grey30")





ggsave(args[[2]], gg, width=21, height=10, units="cm")
