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
keep_subsets <- c(
    #"AllTandemRepeatsandHomopolymers_slop5",
    "AllTandemRepeats_lt51bp_slop5",
    "AllTandemRepeats_51to200bp_slop5",
    "AllTandemRepeats_201to10000bp_slop5",
    "AllTandemRepeats_gt10000bp_slop5"
)

dt <- dt[Subset %in% keep_subsets]
dt[, method := paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(dt$Subset,
                    #"AllTandemRepeatsandHomopolymers_slop5" = "Homopolymers",
                    "AllTandemRepeats_lt51bp_slop5" = "<51",
                    "AllTandemRepeats_51to200bp_slop5" = "51-200",
                    "AllTandemRepeats_201to10000bp_slop5" = "201-10000",
                    "AllTandemRepeats_gt10000bp_slop5" = ">10000"
)

lvls <- c("<51", "51-200", "201-10000",">10000")
dt$Subset <- factor(dt$Subset, levels = lvls)
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

make_heat <- function(df) {
    df %>%
        group_by(cell_line, platform, Subset, tool) %>%
        summarise(
            METRIC.F1_Score = mean(METRIC.F1_Score, na.rm = TRUE),
            QUERY.TOTAL = sum(QUERY.TOTAL, na.rm = TRUE),   # or max(), see note below
            .groups = "drop"
        )
}

snp_h <- make_heat(snp)
idl_h <- make_heat(idl)

layers <- list(
    geom_tile(col = "white", linewidth = 0.1),
    geom_text(aes(label = QUERY.TOTAL), size = 2.5),
    scale_fill_gradientn(
        "F1 score",
        colors = c("ivory", "gold", "red", "navy"),
        na.value = "lightgrey",
        limits = c(0, 1),
        n.breaks = 2
    ),
    theme_minimal(),
    facet_grid2(cell_line ~ platform, scales = "free_y"),
    scale_y_reordered(sep = "___"),
    theme(
        plot.margin = margin(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.text = element_markdown(),
        plot.tag = element_text(size = 9, face = "bold")
    ),
    labs(x="Tandem Repeats length", y="Variant Caller")
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


ggsave(args[[2]], gg, width=32, height=30, units="cm")

