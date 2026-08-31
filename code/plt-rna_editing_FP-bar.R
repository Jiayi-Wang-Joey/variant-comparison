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

dt[, coverage := as.integer(coverage)]
dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
dt[tech == "Baylor-IsoSeq", tech := "IsoSeq"]
dt[, cell_line := sub("^([^-]+)-.*$", "\\1", sample)]
dt[cell_line == "HG002a", cell_line := "HG002"]
dt[, tech := ifelse(grepl("MasSeq|IsoSeq", tech),
                     paste0("<span style='color:#54278f;'>", tech, "</span>"),
                     paste0("<span style='color:#d95f0e;'>", tech, "</span>"))]

d5 <- dt[coverage == 5 & cell_line %in% c("HG002", "HG004", "HG005") & tool != "GATK"]

agg <- d5[, .(FP = sum(FP), FP_in_REDIportal = sum(FP_in_REDIportal_strict)),
          by = .(tool, cell_line, tech)]
agg[, FP_other := FP - FP_in_REDIportal]
agg[, pct := round(100 * FP_in_REDIportal / FP, 1)]
agg[, cell_line := factor(cell_line, levels = c("HG002", "HG004", "HG005"))]

long <- melt(agg, id.vars = c("tool", "cell_line", "tech", "pct", "FP"),
             measure.vars = c("FP_in_REDIportal", "FP_other"),
             variable.name = "category", value.name = "count")
long[, category := factor(category, levels = c("FP_other", "FP_in_REDIportal"),
                           labels = c("Other FP", "RNA editing (REDIportal)"))]

cat_cols <- c("Other FP" = "grey70", "RNA editing (REDIportal)" = "#2166AC")

gg <- ggplot(long, aes(tool, count, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(data = agg, aes(x = tool, y = FP, label = paste0(pct, "%")),
              inherit.aes = FALSE, vjust = -0.3, size = 2) +
    facet_grid2(cell_line ~ tech, scales = "free") +
    theme_classic() +
    scale_fill_manual(values = cat_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = "Variant Caller", y = "Number of FP SNVs", fill = "Category") +
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
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    )

ggsave(args[[2]], gg, width = 30, height = 18, units = "cm")
