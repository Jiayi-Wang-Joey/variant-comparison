suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(ggtext)
})

files <- args[[1]]
files <- files[!grepl("/spectrum/", files)]
dt <- rbindlist(lapply(files, fread))
dt <- dt[Type == "SNP"]

dt[, coverage := factor(as.integer(coverage), levels = sort(unique(as.integer(coverage))))]
dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
dt[tech == "Baylor-IsoSeq", tech := "IsoSeq"]
dt[, cell_line := sub("^([^-]+)-.*$", "\\1", sample)]
dt[cell_line == "HG002a", cell_line := "HG002"]
dt[, tech := ifelse(grepl("MasSeq|IsoSeq", tech),
                     paste0("<span style='color:#54278f;'>", tech, "</span>"),
                     paste0("<span style='color:#d95f0e;'>", tech, "</span>"))]

d <- dt[cell_line %in% c("HG002", "HG004", "HG005") & tool != "GATK"]

agg <- d[, .(FP = sum(FP), FP_in_REDIportal = sum(FP_in_REDIportal_strict)),
         by = .(tool, cell_line, tech, coverage)]
agg[, pct := 100 * FP_in_REDIportal / FP]
agg[, cell_line := factor(cell_line, levels = c("HG002", "HG004", "HG005"))]

tool_cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F"
)

gg <- ggplot(agg, aes(coverage, pct, color = tool, group = tool)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    theme_classic() +
    scale_color_manual(values = tool_cols) +
    labs(x = "Coverage Cutoff (DP >= n)",
         y = "% FP matching REDIportal (position + substitution)",
         color = "Variant Caller") +
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
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    )

ggsave(args[[2]], gg, width = 30, height = 18, units = "cm")
