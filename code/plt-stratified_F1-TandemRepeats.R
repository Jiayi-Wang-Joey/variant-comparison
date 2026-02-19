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

cols <- c("Filter", "METRIC.F1_Score", "coverage", "Subset", "Subtype",
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
keep_subsets <- c(
    #"AllTandemRepeatsandHomopolymers_slop5",
    "AllTandemRepeats_le50bp_slop5",
    "AllTandemRepeats_51to200bp_slop5",
    "AllTandemRepeats_201to10000bp_slop5",
    "AllTandemRepeats_gt10000bp_slop5"
    # "SimpleRepeat_diTR_10to49_slop5",
    # "SimpleRepeat_diTR_50to149_slop5",
    # "SimpleRepeat_triTR_14to49_slop5",	
    # "SimpleRepeat_triTR_50to149_slop5",	
    # "SimpleRepeat_quadTR_19to49_slop5", 
    # "SimpleRepeat_quadTR_50to149_slop5"
    
)

dt <- dt[Subset %in% keep_subsets]
write.table(dt, "data/results/STR.csv")
dt[, method := paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(dt$Subset,
                    "AllTandemRepeats_le50bp_slop5" = "<50",
                    "AllTandemRepeats_51to200bp_slop5" = "51-200",
                    "AllTandemRepeats_201to10000bp_slop5" = "201-10000",
                    "AllTandemRepeats_gt10000bp_slop5" = ">10000"
)
lvls <- c("<50","51-200","201-10000",">10000")
# lvls <- c("diTR_10-49", "diTR_50-149", "triTR_14-49",
#           "triTR_50-149","quadTR_19-49","quadTR_50-149")
dt$Subset <- factor(dt$Subset, levels = lvls)
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[aligner=="minimap2"]
nk <- length(unique(dt$tool))
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk),
                 unique(dt$tool))

snp <- dt[Type=="SNP"]
idl <- dt[Type=="INDEL"]
aes <- list(geom_point(alpha=0.6),
            geom_line(alpha = 0.6),
            facet_grid(cell_line ~ platform),
            theme_minimal(),
            labs(x = "Tandem Repeat length", y = "F1 Score", color = "Variant Caller"),
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
                  strip.text =  element_markdown(),
                  axis.line = element_line(color = "black", linewidth = 0.3),
                  panel.spacing = unit(0, "lines"),
                  panel.spacing.x = unit(0, "lines"),
                  panel.spacing.y = unit(0, "lines"),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))

p1 <- ggplot(snp, aes(Subset, METRIC.F1_Score, 
                      col=tool,
                      group = method)) + aes + ggtitle("SNP") 


p2 <- ggplot(idl, aes(Subset, METRIC.F1_Score, 
                      col=tool,
                      group = method)) + aes + ggtitle("INDEL") + 
    theme(legend.position = "none")

gg <- p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
#write.table(dt, "data/results/STR.csv")
ggsave(args[[2]], gg, width=28, height=25, units="cm")
