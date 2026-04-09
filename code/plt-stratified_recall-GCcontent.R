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

cols <- c("Filter", "METRIC.Recall", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.Recall)]
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
nk <- length(unique(dt$tool))
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk),
                 unique(dt$tool))

snp <- dt[Type=="SNP"]
idl <- dt[Type=="INDEL"]

aes <- list(geom_point(alpha=0.6),
            geom_line(alpha = 0.6),
            facet_grid(cell_line ~ platform),
            theme_minimal(),
            labs(x = "GC Content (%)", y = "Recall", color = "Variant Caller"),
            scale_color_manual(values = cols),
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
                  axis.text.x = element_text(angle = 45, size = 7,
                                             hjust = 1, vjust = 1)))

p1 <- ggplot(snp, aes(Subset, METRIC.Recall, 
                      col=tool,
                      group = method)) + aes + ggtitle("SNP") 


p2 <- ggplot(idl, aes(Subset, METRIC.Recall, 
                      col=tool,
                      group = method)) + aes + ggtitle("INDEL") + 
    theme(legend.position = "none")

gg <- p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold")) 


ggsave(args[[2]], gg, width=32, height=30, units="cm")