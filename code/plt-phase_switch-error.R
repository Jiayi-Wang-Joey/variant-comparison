# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/switch_error/summary/",
# pattern = "\\.tsv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/phase_switch-error.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
    library(scales)
    
})


dt <- rbindlist(lapply(args[[1]], \(f) {
    dt <- fread(f)
    
    keep <- c("sample","phaser","caller","type","n_actual", "n_phased",
              "all_switches","all_assessed_pairs","all_switch_rate")
    dt <- dt[, ..keep]
    
    dt_sum <- dt[all_assessed_pairs > 0,
                 .(
                     all_switches = sum(all_switches),
                     all_assessed_pairs = sum(all_assessed_pairs),
                     mean_switch_rate = sum(all_switches) / sum(all_assessed_pairs),
                     n_chr_informative = .N
                 ),
                 by = .(sample, phaser, caller, type, n_actual, n_phased)
    ]
    
    dt_sum
}), fill = TRUE)

dt[, sample := gsub("-Baylor", "", sample)]

dt <- dt[grepl("^HG002a-|^HG005-", sample) | sample == "HG002-dRNA004"]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, platform := factor(platform, levels = sort(unique(platform)))]

dt[, cell_line := sapply(strsplit(sample, "-"), head, 1)]

dt[cell_line == "HG002a", cell_line := "HG002"]
dt[, cell_line := factor(cell_line)]

dt <- dt[caller != "GATK"]

dt[, method := paste0(caller, phaser)]
dt[phaser == "longcallR", phaser := "longcallR*"]
dt[phaser == "isoLASER",  phaser := "isoLASER*"]

snp <- dt[type == "SNP" & caller != "isoLASER" & phaser != "longcallR*"]

all <- dt[type == "all" & !grepl("longcallR|isoLASER|NanoTS", caller)]

# longcallR*/isoLASER*: repeated in every caller row
callers_in_plot <- sort(unique(snp$caller))

snp_star <- dt[type == "SNP" & phaser %in% c("longcallR*", "isoLASER*"),
               .(sample, phaser, type, platform, all_assessed_pairs,
                 mean_switch_rate, cell_line)]

star <- rbindlist(lapply(callers_in_plot, \(cl) {
    copy(snp_star)[, caller := cl]
}))

# isoLASER*: PacBio only
callers_in_plot_all <- sort(unique(all$caller))

all_star <- dt[type == "all" & phaser == "isoLASER*" & grepl("MasSeq|IsoSeq", sample),
               .(sample, phaser, type, platform, all_assessed_pairs,
                 mean_switch_rate, cell_line)]

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
    facet_grid2(caller ~ platform, scales = "free", drop = FALSE),
    theme_classic(),
    scale_x_continuous(labels = label_scientific(digits = 2)),
    labs(y = "1 - Switch Error Rate", x = "Number of Assessed Pairs",
         shape = "Cell Line", color = "Phaser"),
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

p1 <- ggplot(snp, aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
                      color = phaser, shape = cell_line)) +
    geom_point(
        data = star,
        aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
            color = phaser, shape = cell_line),
        size = 3, alpha = 0.8,
        inherit.aes = FALSE
    ) +
    aes_common +
    ggtitle("SNV")

p2 <- ggplot(all, aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
                      color = phaser, shape = cell_line)) +
    geom_point(
        data = star_all,
        aes(x = all_assessed_pairs, y = 1 - mean_switch_rate,
            color = phaser, shape = cell_line),
        size = 3, alpha = 0.8,
        inherit.aes = FALSE
    ) +
    aes_common +
    guides(color = "none", shape = "none") +
    ggtitle("ALL")

gg <- (p1 / p2) +
    plot_layout(guides = "collect", widths = c(1, 1), heights = c(2, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width = 22, height = 27, units = "cm")