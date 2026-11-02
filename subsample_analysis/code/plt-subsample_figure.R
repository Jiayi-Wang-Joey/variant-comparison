suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(ggtext)
})

# y-positions for stacked-bar count labels, nudged apart to avoid overlap
add_label_pos <- function(dt, by, min_gap_frac = 0.18) {
    dt <- copy(dt)
    setorderv(dt, c(by, "Category"))
    dt[, ymax := rev(cumsum(rev(Count))), by = by]
    dt[, ymin := ymax - Count]
    dt[, ymid := (ymin + ymax) / 2]
    dt[, total := sum(Count), by = by]
    dt[, label_y := {
        o <- order(ymid)
        y <- ymid[o]
        mg <- total[1] * min_gap_frac
        n <- length(y)
        if (n > 1) {
            for (i in 2:n) if (y[i] - y[i - 1] < mg) y[i] <- y[i - 1] + mg
        }
        out <- numeric(n)
        out[o] <- y
        out
    }, by = by]
    dt
}

args <- commandArgs(trailingOnly = TRUE)
bases_files <- strsplit(args[1], ";")[[1]]
happy_files <- strsplit(args[2], ";")[[1]]
out_pdf <- args[3]

SOURCE_SAMPLE <- "HG004-MasSeq"
REF_SAMPLE <- "HG004-dRNA004"
SUBSAMPLE_SEEDS <- 1:5
SUBSAMPLE_SAMPLES <- paste0("HG004-MasSeq-subsampled-seed", SUBSAMPLE_SEEDS)


sample_levels <- c(SOURCE_SAMPLE, SUBSAMPLE_SAMPLES, REF_SAMPLE)
sample_label_of <- c(
    setNames("MasSeq (origin)", SOURCE_SAMPLE),
    setNames(paste0("MasSeq (seed ", SUBSAMPLE_SEEDS, ")"), SUBSAMPLE_SAMPLES),
    setNames("dRNA004", REF_SAMPLE)
)
seed_cols <- colorRampPalette(c("#cbc9e2", "#54278f"))(length(SUBSAMPLE_SEEDS))
sample_cols <- c(
    setNames("#238b45", SOURCE_SAMPLE),
    setNames(seed_cols, SUBSAMPLE_SAMPLES),
    setNames("#d95f0e", REF_SAMPLE)
)

sample_scale <- scale_fill_manual(
    values = sample_cols, aesthetics = c("fill", "colour"), name = "",
    breaks = sample_levels, labels = sample_label_of[sample_levels]
)

# ---- panel a -----
bases <- rbindlist(lapply(bases_files, fread, header = TRUE))
bases <- bases[metric == "callable_bases"]
bases[, sample := factor(sample, levels = sample_levels)]
bases[, coverage := factor(min_coverage, levels = sort(unique(min_coverage)))]

p_a <- ggplot(bases, aes(coverage, value, fill = sample, group = sample)) +
    geom_col(position = position_dodge2(width = 0.8, padding = 0.1)) +
    sample_scale +
    labs(x = "Coverage Cutoff (DP >= n)", y = "Callable bases") +
    theme_classic() +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom"
    ) 

# ---- panel b ---------
happy <- rbindlist(lapply(happy_files, fread, header = TRUE))
happy[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]

happy_pass <- happy[Filter == "PASS" & !is.na(METRIC.F1_Score)]
happy_pass[, sample := factor(sample, levels = sample_levels)]
happy_pass[, coverage := factor(coverage, levels = sort(unique(coverage)))]

p_b <- ggplot(happy_pass, aes(coverage, METRIC.F1_Score, color = sample, group = sample)) +
    geom_point(size = 1.5, alpha = 0.9) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    sample_scale +
    guides(colour = "none", fill = "none") +
    facet_grid2(. ~ Type, scales = "free_y", independent = "y") +
    labs(x = "Coverage Cutoff (DP >= n)", y = "F1 Score") +
    theme_classic() +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5)
    )

# ---- panel c ------------
PANEL_C_COVERAGE <- c(5, 100)
snapshot <- happy[coverage %in% PANEL_C_COVERAGE & Filter == "PASS"]
snapshot[, coverage_row := factor(
    paste0("DP>=", coverage), levels = paste0("DP>=", PANEL_C_COVERAGE)
)]

origin_dt <- snapshot[sample %in% c(SOURCE_SAMPLE, REF_SAMPLE)]
origin_dt[, facet_col := "origin"]
origin_dt[, sample_label := ifelse(sample == SOURCE_SAMPLE, "MasSeq", "dRNA004")]

subsampled_dt <- snapshot[sample %in% SUBSAMPLE_SAMPLES]
subsampled_dt[, facet_col := "subsampled"]
subsampled_dt[, sample_label := sub("HG004-MasSeq-subsampled-seed", "seed", sample)]

panel_c_dt <- rbindlist(list(origin_dt, subsampled_dt))
panel_c_dt[, facet_col := factor(facet_col, levels = c("origin", "subsampled"))]
panel_c_dt[, sample_label := factor(
    sample_label, levels = c("MasSeq", paste0("seed", SUBSAMPLE_SEEDS), "dRNA004")
)]

td <- melt(
    panel_c_dt,
    id.vars = c("Type", "facet_col", "sample_label", "coverage_row"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP", "METRIC.F1_Score"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]

.panel_c <- function(dt, title) {
    counts <- dt[Category %in% c("TP", "FN", "FP")]
    counts_lab <- add_label_pos(counts, by = c("facet_col", "sample_label", "coverage_row"))
    stack_top <- counts_lab[, .(top_y = max(label_y), total = total[1]),
                             by = .(facet_col, sample_label, coverage_row)]
    f1 <- merge(
        dt[Category == "F1_Score", .(facet_col, sample_label, coverage_row, F1 = Count)],
        stack_top,
        by = c("facet_col", "sample_label", "coverage_row")
    )
    f1[, label_y := top_y + total * 0.18]
    ggplot() +
        geom_col(data = counts, aes(sample_label, Count, fill = Category)) +
        geom_text(
            data = counts_lab,
            aes(sample_label, label_y, label = Count),
            size = 2.4
        ) +
        geom_text(
            data = f1,
            aes(sample_label, label_y, label = sprintf("%.2f", F1)),
            size = 2.8, fontface = "bold"
        ) +
        facet_grid2(coverage_row ~ facet_col, scales = "free", space = "free_x") +
        scale_fill_brewer(palette = "Paired") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
        labs(x = "", y = "Total") +
        theme_classic() +
        theme(
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
        ) +
        ggtitle(title)
}

p_c1 <- .panel_c(td[Type == "SNV"], "SNV")
p_c2 <- .panel_c(td[Type == "INDEL"], "INDEL")


p_ab <- (p_a / p_b) +
    plot_layout(heights = c(1, 1.5))

p_cd <- (p_c1 + p_c2) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.box.spacing = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0)
    )

gg <- wrap_plots(p_ab, p_cd, ncol = 1, heights = c(1.6, 1.4)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(out_pdf, gg, width = 22, height = 28, units = "cm")
