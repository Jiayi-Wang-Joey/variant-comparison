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

cols_keep <- c("Filter", "METRIC.F1_Score", "coverage", "Subset", "Subtype",
               "sample", "aligner", "tool", "bamtype", "Type", "TRUTH.TOTAL",
               "TRUTH.TP", "QUERY.TOTAL")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols_keep] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter == "PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR", tool) & Type == "INDEL")]
dt <- dt[coverage == 5 & aligner == "minimap2"]

dt[, cell_line := sapply(strsplit(sample, "-"), head, 1)]
dt[cell_line == "HG002a", cell_line := "HG002"]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, method := paste(aligner, bamtype, tool, sep = ".")]

saveRDS(dt, "data/results/stratified.rds")

cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F"
)

theme_ctx <- theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_markdown(size = 11),
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.spacing = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines"),
    panel.spacing.y = unit(0, "lines"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, size = 7, hjust = 1, vjust = 1),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    legend.title = element_text(size = 11)
)

.panel <- \(d, x_lab, title) {
    d[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]
    lb <- d[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = c("Subset", "platform", "Type")]

    ggplot(d, aes(Subset, METRIC.F1_Score, color = tool, group = method)) +
        geom_point(alpha = 0.6) +
        geom_line(alpha = 0.6) +
        geom_text(data = lb,
            aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
            inherit.aes = FALSE,
            angle = 30, hjust = 0.5, vjust = 1.3,
            size = 2, color = "grey30") +
        coord_cartesian(ylim = c(0, 1.15), clip = "off") +
        facet_grid(Type ~ platform) +
        theme_classic() +
        labs(x = x_lab, y = "F1 Score", color = "Variant Caller") +
        scale_color_manual(values = cols) +
        theme_ctx +
        ggtitle(title)
}

tr_subsets <- c(
    "SimpleRepeat_diTR_10to49_slop5",
    "SimpleRepeat_diTR_50to149_slop5",
    "SimpleRepeat_triTR_14to49_slop5",
    "SimpleRepeat_triTR_50to149_slop5",
    "SimpleRepeat_quadTR_19to49_slop5",
    "SimpleRepeat_quadTR_50to149_slop5"
)
gc_lvls <- c("<15", "15-20", "20-25", "25-30", "30-55", "55-60",
             "60-65", "65-70", "70-75", "75-80", "80-85", ">85")

.page <- \(cl) {
    d <- dt[cell_line == cl]

    ## ---- Homopolymer ----
    homo <- d[grepl("homopolymer", Subset) & Subtype == "*"]
    homo <- homo[!grepl("AT|GC", Subset)]
    homo$Subset <- recode(homo$Subset,
        "SimpleRepeat_homopolymer_4to6_slop5"  = "4–6",
        "SimpleRepeat_homopolymer_7to11_slop5" = "7–11",
        "SimpleRepeat_homopolymer_ge12_slop5"  = "≥12",
        "SimpleRepeat_homopolymer_ge21_slop5"  = "≥21"
    )
    homo$Subset <- factor(homo$Subset, levels = c("4–6", "7–11", "≥12", "≥21"))
    panel_homo <- .panel(homo, "Homopolymer length", paste0(cl, " — Homopolymer"))

    ## ---- Tandem Repeat ----
    tr <- d[Subset %in% tr_subsets & Subtype == "*"]
    tr$Subset <- recode(tr$Subset,
        "SimpleRepeat_diTR_10to49_slop5"    = "diTR_10-49",
        "SimpleRepeat_diTR_50to149_slop5"   = "diTR_50-149",
        "SimpleRepeat_triTR_14to49_slop5"   = "triTR_14-49",
        "SimpleRepeat_triTR_50to149_slop5"  = "triTR_50-149",
        "SimpleRepeat_quadTR_19to49_slop5"  = "quadTR_19-49",
        "SimpleRepeat_quadTR_50to149_slop5" = "quadTR_50-149"
    )
    tr$Subset <- factor(tr$Subset, levels = c(
        "diTR_10-49", "diTR_50-149", "triTR_14-49",
        "triTR_50-149", "quadTR_19-49", "quadTR_50-149"
    ))
    panel_tr <- .panel(tr, "Tandem Repeat Subclass", "Tandem Repeat")

    ## ---- GC content ----
    gc <- d[grepl("^gc[0-9]", Subset) & Subtype == "*"]
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
    gc$Subset <- factor(gc$Subset, levels = gc_lvls)
    panel_gc <- .panel(gc, "GC Content (%)", "GC content")

    panel_homo / panel_tr / panel_gc +
        plot_layout(heights = c(1, 1, 1), guides = "collect") +
        plot_annotation(tag_levels = list(c("b", "c", "d"))) &
        theme(plot.tag = element_text(face = "bold"), legend.position = "bottom")
}

pdf(args[[2]], width = 22 / 2.54, height = 30 / 2.54)
for (cl in c("HG004", "HG002", "HG005")) {
    print(.page(cl))
}
dev.off()
