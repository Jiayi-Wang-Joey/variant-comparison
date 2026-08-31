suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
    library(scales)
})

dt <- rbindlist(lapply(args[[1]], fread), fill = TRUE)

dt[, sample := gsub("-Baylor", "", sample)]

# Drop stale pre-"HG002a" naming leftovers (e.g. "HG002-MasSeq",
# "HG002-Baylor-IsoSeq") from old runs — only HG002a-* and the
# intentionally-plain HG002-dRNA004 are part of the current sample set.
dt <- dt[grepl("^HG002a-|^HG005-", sample) | sample == "HG002-dRNA004"]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := sapply(strsplit(sample, "-"), head, 1)]

# Most assays are labelled "HG002a", but dRNA004 is already plain "HG002" —
# both refer to the same cell line, so unify the label.
dt[cell_line == "HG002a", cell_line := "HG002"]
dt[, cell_line := factor(cell_line)]

dt[phaser == "longcallR", phaser := "longcallR*"]
dt[phaser == "isoLASER",  phaser := "isoLASER*"]

# Main SNP data: only standard phasers (have caller facet rows)
snp <- dt[type == "SNP" & caller != "isoLASER" & phaser != "longcallR*"]

# ALL data: exclude longcallR/isoLASER callers
all <- dt[type == "all" & !grepl("longcallR|isoLASER", caller)]

# isoLASER* and longcallR* SNP data — replicated across all caller facet rows
callers_in_plot <- sort(unique(snp$caller))

snp_star <- dt[type == "SNP" & phaser %in% c("longcallR*", "isoLASER*")]

star <- rbindlist(lapply(callers_in_plot, \(cl) {
    copy(snp_star)[, caller := cl]
}))

# isoLASER* is PacBio-only — replicate across ALL-panel caller facet rows,
# but only for PacBio platforms (IsoSeq/MasSeq)
callers_in_plot_all <- sort(unique(all$caller))
all_star <- dt[type == "all" & phaser == "isoLASER*" & grepl("MasSeq|IsoSeq", sample)]

star_all <- rbindlist(lapply(callers_in_plot_all, \(cl) {
    copy(all_star)[, caller := cl]
}))

cols <- c(
    "WhatsHap"   = "#1b9e77",
    "HapCUT2"    = "#d95f02",
    "longphase"  = "#377eb8",
    "HiPhase"    = "#e7298a",
    "longcallR*" = "grey60",
    "isoLASER*"  = "#984ea3"
)

aes_common <- list(
    geom_point(size = 3, alpha = 0.8),
    facet_grid2(caller ~ platform, scales = "free"),
    theme_classic(),
    labs(shape = "Cell Line", color = "Phaser"),
    theme(
        panel.grid.major   = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.8),
        strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text         = element_markdown(size = 11),
        axis.line          = element_line(color = "black", linewidth = 0.3),
        panel.spacing      = unit(0, "lines"),
        panel.spacing.x    = unit(0, "lines"),
        panel.spacing.y    = unit(0, "lines"),
        axis.text.y        = element_text(size = 7),
        axis.title.x       = element_text(size = 11),
        axis.title.y       = element_text(size = 11),
        legend.title       = element_text(size = 11),
        axis.text.x        = element_text(angle = 45, size = 7, hjust = 1, vjust = 1),
        aspect.ratio       = 1
    ),
    scale_color_manual(values = cols)
)

# x = fraction of het. variants phased, y = block N50 (bp)
.mkplot <- \(d, s, show_legend = TRUE) {
    p <- ggplot(d, aes(x = phased_fraction, y = block_n50_bp,
                       color = phaser, shape = cell_line)) +
        geom_point(
            data = s,
            aes(x = phased_fraction, y = block_n50_bp,
                color = phaser, shape = cell_line),
            size = 3, alpha = 0.8,
            inherit.aes = FALSE
        ) +
        aes_common +
        scale_y_continuous(trans = "log10", labels = label_scientific(digits = 2)) +
        labs(x = "Fraction of Het. Variants Phased", y = "Phase Block N50 (bp)")
    if (!show_legend) p <- p + guides(color = "none", shape = "none")
    p
}

n50_snp <- .mkplot(snp, star) + ggtitle("SNV")
n50_all <- .mkplot(all, star_all, show_legend = FALSE) + ggtitle("ALL")

gg <- (n50_snp / n50_all) +
    plot_layout(guides = "collect", heights = c(1.5, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width = 22, height = 22, units = "cm")
