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
    library(tidytext)
})

cols <- c("Filter", "METRIC.F1_Score", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type", "QUERY.TOTAL")
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
make_heat <- function(df) {
    df %>%
        group_by(cell_line, platform, Subset, tool,text_col) %>%
        summarise(
            METRIC.F1_Score = mean(METRIC.F1_Score, na.rm = TRUE),
            QUERY.TOTAL = sum(QUERY.TOTAL, na.rm = TRUE),   
            .groups = "drop"
        )
}
dt[, text_col := fifelse(METRIC.F1_Score >= 0.9, 
                         "white", "black")]
snp <- dt[Type=="SNP"]
idl <- dt[Type=="INDEL"]

snp_h <- make_heat(snp)
idl_h <- make_heat(idl)

layers <- list(
    geom_tile(col = "white", linewidth = 0.1),
    geom_text(aes(label = QUERY.TOTAL,
                  color = text_col), 
              size = 1.5),
    scale_fill_gradientn(
        "F1 score",
        colors = c("ivory", "gold", "red", "navy"),
        na.value = "lightgrey",
        limits = c(0, 1),
        n.breaks = 2
    ),
    scale_color_identity(),
    theme_minimal(),
    facet_grid2(cell_line ~ platform, scales = "free_y"),
    scale_y_reordered(sep = "___"),
    theme(
        plot.margin = margin(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.text = element_markdown(),
        plot.tag = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ),
    labs(x="GC Content (%)", y="Variant Caller")
)
SEP <- "|" 
p1 <- ggplot(
    snp_h,
    aes(
        Subset,
        reorder_within(tool, METRIC.F1_Score, cell_line, sep = SEP, desc = TRUE),
        fill = METRIC.F1_Score
    )
) +
    layers +
    scale_y_discrete(labels = function(x) sub("\\|.*$", "", x))

p2 <- ggplot(
    idl_h,
    aes(
        Subset,
        reorder_within(tool, METRIC.F1_Score, cell_line, sep = SEP, desc = TRUE),
        fill = METRIC.F1_Score
    )
) +
    layers +
    scale_y_discrete(labels = function(x) sub("\\|.*$", "", x))

gg <- p1 + p2 +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))


ggsave(args[[2]], gg, width=50, height=25, units="cm")
write.table(dt, "data/results/GC_content.csv")