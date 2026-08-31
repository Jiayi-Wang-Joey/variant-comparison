suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
    library(scales)
})

files <- args[[1]]
se_files <- files[grepl("switch_error", files)]
pb_files <- files[grepl("phase_blocks", files)]

# --- accuracy data (switch error / assessed pairs): per-chromosome rows,
# aggregate to one row per sample/caller/phaser/type
dt_se <- rbindlist(lapply(se_files, \(f) {
    d <- fread(f)
    keep <- c("sample", "phaser", "caller", "type", "n_actual", "n_phased",
              "all_switches", "all_assessed_pairs", "all_switch_rate")
    d <- d[, ..keep]
    d[all_assessed_pairs > 0,
      .(all_switches = sum(all_switches),
        all_assessed_pairs = sum(all_assessed_pairs),
        mean_switch_rate = sum(all_switches) / sum(all_assessed_pairs),
        n_chr_informative = .N),
      by = .(sample, phaser, caller, type, n_actual, n_phased)]
}), fill = TRUE)

# --- coverage data (phase block N50 / fraction phased): already one row per combo
dt_pb <- rbindlist(lapply(pb_files, fread), fill = TRUE)

# --- shared munging (sample/platform/cell_line/phaser labels), applied to both
.munge <- \(dt) {
    dt[, sample := gsub("-Baylor", "", sample)]

    # Drop stale pre-"HG002a" naming leftovers (e.g. "HG002-MasSeq",
    # "HG002-Baylor-IsoSeq") from old runs — only HG002a-* and the
    # intentionally-plain HG002-dRNA004 are part of the current sample set.
    dt <- dt[grepl("^HG002a-|^HG005-", sample) | sample == "HG002-dRNA004"]

    dt[, platform := sapply(strsplit(sample, "-"), tail, 1)]
    dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                            paste0("<span style='color:#54278f;'>", platform, "</span>"),
                            paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]

    dt[, cell_line := sapply(strsplit(sample, "-"), head, 1)]
    # Most assays are labelled "HG002a", but dRNA004 is already plain "HG002" —
    # both refer to the same cell line, so unify the label.
    dt[cell_line == "HG002a", cell_line := "HG002"]
    dt[, cell_line := factor(cell_line)]

    dt <- dt[caller != "GATK"]
    dt[phaser == "longcallR", phaser := "longcallR*"]
    dt[phaser == "isoLASER",  phaser := "isoLASER*"]
    dt
}

dt_se <- .munge(dt_se)
dt_pb <- .munge(dt_pb)

# Shared platform factor levels across both datasets, so all four panels
# (accuracy SNV/ALL, coverage SNV/ALL) render the same facet columns and stay
# the same width when stacked.
platform_levels <- sort(unique(c(dt_se$platform, dt_pb$platform)))
dt_se[, platform := factor(platform, levels = platform_levels)]
dt_pb[, platform := factor(platform, levels = platform_levels)]

cols <- c(
    "WhatsHap"   = "#1b9e77",
    "HapCUT2"    = "#d95f02",
    "longphase"  = "#377eb8",
    "HiPhase"    = "#e7298a",
    "longcallR*" = "grey60",
    "isoLASER*"  = "#984ea3"
)

theme_common <- list(
    geom_point(size = 3, alpha = 0.8),
    facet_grid2(caller ~ platform, scales = "free", drop = FALSE),
    theme_classic(),
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
        aspect.ratio       = 1,
        plot.title         = element_text(size = 11, margin = margin(b = 2)),
        plot.margin        = margin(t = 2, r = 2, b = 2, l = 2)
    ),
    scale_color_manual(values = cols)
)

# ---- accuracy panels: switch error rate vs number of assessed pairs ----
snp_se <- dt_se[type == "SNP" & caller != "isoLASER" & phaser != "longcallR*"]
all_se <- dt_se[type == "all" & !grepl("longcallR|isoLASER", caller)]

callers_se <- sort(unique(snp_se$caller))
snp_star_se <- dt_se[type == "SNP" & phaser %in% c("longcallR*", "isoLASER*"),
                     .(sample, phaser, type, platform, all_assessed_pairs,
                       mean_switch_rate, cell_line)]
star_se <- rbindlist(lapply(callers_se, \(cl) copy(snp_star_se)[, caller := cl]))

callers_all_se <- sort(unique(all_se$caller))
all_star_se <- dt_se[type == "all" & phaser == "isoLASER*" & grepl("MasSeq|IsoSeq", sample),
                     .(sample, phaser, type, platform, all_assessed_pairs,
                       mean_switch_rate, cell_line)]
star_all_se <- rbindlist(lapply(callers_all_se, \(cl) copy(all_star_se)[, caller := cl]))

p1 <- ggplot(snp_se, aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
                         color = phaser, shape = cell_line)) +
    geom_point(
        data = star_se,
        aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
            color = phaser, shape = cell_line),
        size = 3, alpha = 0.8, inherit.aes = FALSE
    ) +
    theme_common +
    scale_x_continuous(labels = label_scientific(digits = 2)) +
    labs(y = "1 - Switch Error Rate", x = "Number of Assessed Pairs",
         shape = "Cell Line", color = "Phaser") +
    ggtitle("SNV")

p2 <- ggplot(all_se, aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
                         color = phaser, shape = cell_line)) +
    geom_point(
        data = star_all_se,
        aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
            color = phaser, shape = cell_line),
        size = 3, alpha = 0.8, inherit.aes = FALSE
    ) +
    theme_common +
    scale_x_continuous(labels = label_scientific(digits = 2)) +
    labs(y = "1 - Switch Error Rate", x = "Number of Assessed Pairs",
         shape = "Cell Line", color = "Phaser") +
    guides(color = "none", shape = "none") +
    ggtitle("ALL")

# ---- coverage panels: phase block N50 vs fraction of het. variants phased ----
snp_pb <- dt_pb[type == "SNP" & caller != "isoLASER" & phaser != "longcallR*"]
all_pb <- dt_pb[type == "all" & !grepl("longcallR|isoLASER", caller)]

callers_pb <- sort(unique(snp_pb$caller))
snp_star_pb <- dt_pb[type == "SNP" & phaser %in% c("longcallR*", "isoLASER*")]
star_pb <- rbindlist(lapply(callers_pb, \(cl) copy(snp_star_pb)[, caller := cl]))

callers_all_pb <- sort(unique(all_pb$caller))
all_star_pb <- dt_pb[type == "all" & phaser == "isoLASER*" & grepl("MasSeq|IsoSeq", sample)]
star_all_pb <- rbindlist(lapply(callers_all_pb, \(cl) copy(all_star_pb)[, caller := cl]))

.mkplot_pb <- \(d, s, show_legend = TRUE) {
    p <- ggplot(d, aes(x = phased_fraction, y = block_n50_bp,
                       color = phaser, shape = cell_line)) +
        geom_point(
            data = s,
            aes(x = phased_fraction, y = block_n50_bp,
                color = phaser, shape = cell_line),
            size = 3, alpha = 0.8, inherit.aes = FALSE
        ) +
        theme_common +
        scale_y_continuous(trans = "log10", labels = label_scientific(digits = 2)) +
        labs(x = "Fraction of Het. Variants Phased", y = "Phase Block N50 (bp)",
             shape = "Cell Line", color = "Phaser")
    if (!show_legend) p <- p + guides(color = "none", shape = "none")
    p
}

n50_snp <- .mkplot_pb(snp_pb, star_pb) + ggtitle("SNV")
n50_all <- .mkplot_pb(all_pb, star_all_pb, show_legend = FALSE) + ggtitle("ALL")

# ---- combine as a flat 2x2 grid: coverage (N50) column on the left
# (a=SNV, b=ALL), accuracy (switch error) column on the right (c=SNV,
# d=ALL). Using wrap_plots(..., byrow=FALSE) instead of nested "/"/"|"
# so `heights` unambiguously applies to the two grid rows — SNV has 3
# caller rows vs ALL's 2, and with aspect.ratio=1 each panel cell is
# forced square, so the row heights must match that 3:2 ratio or the
# shorter (ALL) panels are padded with whitespace to match a shared height.
gg <- wrap_plots(list(n50_snp, n50_all, p1, p2), ncol = 2, byrow = FALSE, heights = c(3, 2)) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width = 38, height = 19, units = "cm")
