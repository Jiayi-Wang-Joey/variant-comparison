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
keep_subsets <- c(
    "SimpleRepeat_diTR_10to49_slop5",
    "SimpleRepeat_diTR_50to149_slop5",
    "SimpleRepeat_diTR_ge150_slop5",
    "SimpleRepeat_triTR_14to49_slop5",	
    "SimpleRepeat_triTR_50to149_slop5",	
    "SimpleRepeat_triTR_ge150_slop5",
    "SimpleRepeat_quadTR_19to49_slop5", 
    "SimpleRepeat_quadTR_50to149_slop5",
    "SimpleRepeat_quadTR_ge150_slop5"
)
dt <- dt[Subset %in% keep_subsets]
dt$Subset <- recode(dt$Subset,
                    "SimpleRepeat_diTR_10to49_slop5"="diTR_10-49",
                    "SimpleRepeat_diTR_50to149_slop5"="diTR_50-149",
                    "SimpleRepeat_diTR_ge150_slop5"="diTR_>150",
                    "SimpleRepeat_triTR_14to49_slop5"="triTR_14-49",
                    "SimpleRepeat_triTR_50to149_slop5"="triTR_50-149",
                    "SimpleRepeat_triTR_ge150_slop5"="triTR_>150",
                    "SimpleRepeat_quadTR_19to49_slop5"="quadTR_19-49",
                    "SimpleRepeat_quadTR_50to149_slop5"="quadTR_50-149",
                    "SimpleRepeat_quadTR_ge150_slop5"="quadTR_>150"
)
write.table(dt, "data/results/STR.csv")
dt[, method := paste(aligner, bamtype, tool, sep = ".")]

lvls <- c(
    "diTR_10-49", "diTR_50-149","diTR_>150",
    "triTR_14-49", "triTR_50-149","triTR_>150",
    "quadTR_19-49", "quadTR_50-149", "quadTR_>150"
)
dt$Subset <- factor(dt$Subset, levels = lvls)
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt[cell_line == "HG002a", cell_line := "HG002"]
dt <- dt[aligner=="minimap2"]
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F"
)
dt <- dt[!(Type == "INDEL" & grepl("dRNA|cDNA", platform) & tool == "isoLASER")]
# HG002a is relabeled to HG002 above; average the now-duplicate HG002/HG002a
# points per (Subset, platform, tool) instead of plotting both separately
dt <- dt[, .(METRIC.F1_Score = mean(METRIC.F1_Score, na.rm = TRUE),
             TRUTH.TOTAL = mean(TRUTH.TOTAL, na.rm = TRUE)),
         by = .(Subset, platform, cell_line, tool, aligner, bamtype, Type)]
dt[, method := paste(aligner, bamtype, tool, sep = ".")]
snp <- dt[Type=="SNP"]
idl <- dt[Type=="INDEL"]
lb1 <- snp[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(Subset, cell_line, platform)]
lb2 <- idl[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(Subset, cell_line, platform)]
aes <- list(geom_point(alpha=0.6),
            geom_line(alpha = 0.6),
            facet_grid(cell_line ~ platform),
            theme_minimal(),
            labs(x = "Tandem Repeat Subclass", y = "F1 Score", color = "Variant Caller"),
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
                  aspect.ratio = 1,
                  axis.text.x = element_text(angle = 45, size = 5,
                                             hjust = 1, vjust = 1)))

p1 <- ggplot(snp, aes(Subset, METRIC.F1_Score,
                      col=tool,
                      group = method)) + aes +
    geom_text(data = lb1,
        aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
        inherit.aes = FALSE,
        angle = 30, hjust = 0.5, vjust = 1.3,
        size = 2.2, color = "grey30") +
    ggtitle("SNV")


p2 <- ggplot(idl, aes(Subset, METRIC.F1_Score,
                      col=tool,
                      group = method)) + aes +
    geom_text(data = lb2,
        aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
        inherit.aes = FALSE,
        angle = 30, hjust = 0.5, vjust = 1.3,
        size = 2.2, color = "grey30") +
    ggtitle("INDEL") +
    theme(legend.position = "none")

gg <- p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))
write.table(dt, "data/results/STR.csv")
ggsave(args[[2]], gg, width=26, height=28, units="cm")
