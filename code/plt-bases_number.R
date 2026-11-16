#args <- list(list.files("~/omni/jiayiwang/data/variant-comparison/results/beds/raw/summary/",
#pattern = "\\.tsv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/bases-number.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
    library(RColorBrewer)
})

# a: standard region, b: common region (HG004)
res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res, fill = TRUE)

dt <- dt[metric == "callable_bases"]
dt[, coverage := factor(min_coverage, levels = unique(sort(min_coverage)))]

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1),
    levels = c("IsoSeq", "MasSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt[cell_line=="HG002a", cell_line := "HG002"]

cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER" = "#FDBF6F"
)

aes <- list(
    theme_classic(),
    facet_grid2(cell_line ~ platform, scales = "free"),
    labs(
        x = "Coverage Cutoff (DP >= n)",
        y = "Callable bases"
    ),
    scale_color_manual(values = cols),
    scale_x_discrete(expand = expansion(add = c(0.8, 1))),
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
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    )
)


.bars <- \(d) ggplot(d, aes(coverage, value / 1e6)) +
    geom_bar(stat = "identity") +
    aes +
    labs(y = "Callable bases (Mb)")

pa <- .bars(dt[status == "raw"]) + labs(title = "Standard evaluation region")
pb <- .bars(dt[status == "common"]) + labs(title = "Common region (standard region with isoLASER effective depth)")
gg <- pa / pb + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(args[[2]], gg, width=25, height=34, units="cm")
