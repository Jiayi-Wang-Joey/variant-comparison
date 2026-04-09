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
    
    keep <- c(
        "sample","phaser","caller","type",
        "n_actual","n_phased",
        "all_switches","all_assessed_pairs",
        "blockwise_hamming"
    )
    dt <- dt[, ..keep]
    
    dt_sum <- dt[all_assessed_pairs > 0,
                 .(
                     all_switches = sum(all_switches),
                     all_assessed_pairs = sum(all_assessed_pairs),
                     mean_switch_rate = sum(all_switches) / sum(all_assessed_pairs),
                     
                     blockwise_hamming = sum(blockwise_hamming),
                     mean_hamming_rate = sum(blockwise_hamming) / sum(all_assessed_pairs),
                     
                     n_chr_informative = .N
                 ),
                 by = .(sample, phaser, caller, type, n_actual, n_phased)
    ]
    
    dt_sum
}), fill = TRUE)

dt[, sample:=gsub("-Baylor", "", sample)]

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]

dt[,method:=paste0(caller, phaser)]
dt[,phaser:=factor(phaser)]
snp <- dt[type=="SNP" & phaser != "longcallR"]
all <- dt[type=="all" & !grepl("longcallR", caller)]
snp_longcallR <- dt[type=="SNP" & phaser == "longcallR"]
callers_in_plot <- sort(unique(snp$caller))
star <- snp_longcallR[
    , .(sample, phaser, type, platform, all_assessed_pairs, mean_hamming_rate, cell_line)
][
    , caller := callers_in_plot[1]
]
star <- star[rep(seq_len(.N), each = length(callers_in_plot))]
star[, caller := rep(callers_in_plot, times = nrow(snp_longcallR))]
cols <- c(
    "WhatsHap"  = "#66c2a5",
    "HapCUT2"   = "#fc8d62",
    "longphase" = "#8da0cb",
    "HiPhase"   = "#e78ac3",
    "longcallR" = "grey70"
)


aes <- list(
    geom_point(size=3, alpha = 0.8),
    facet_grid2(caller ~ platform, scales="free"),
    theme_classic(),
    scale_x_continuous(labels = label_scientific(digits = 2)),
    labs(y = "1 - Hamming Rate", x = "Number of Assessed Pairs",
         shape = "Cell Line", color = "Phaser"),
    theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
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
          strip.text =  element_markdown(size=11),
          axis.line = element_line(color = "black", linewidth = 0.3),
          panel.spacing = unit(0, "lines"),
          panel.spacing.x = unit(0, "lines"),
          panel.spacing.y = unit(0, "lines"),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text.x = element_text(angle = 45, size = 7,
                                     hjust = 1, vjust = 1),
          aspect.ratio = 1),
    #scale_color_brewer(palette = "Set2")
    scale_color_manual(values=cols)
)
p1 <- ggplot(snp, aes(x = all_assessed_pairs, y = 1 - mean_hamming_rate, 
                      color = phaser, shape = cell_line)) + 
    geom_point(
        data = star,
        aes(x = all_assessed_pairs,
            y = 1 - mean_hamming_rate,
            shape = cell_line,
            color = "longcallR"),
        size = 3,
        alpha = 0.8,
        inherit.aes = FALSE
    ) +
    aes + ggtitle("SNP")

p2 <- ggplot(all, aes(x = all_assessed_pairs, y = 1 - mean_hamming_rate, 
                      color = phaser, shape = cell_line)) + aes + 
    guides(color = "none",
           shape = "none") + ggtitle("ALL")


# p1 <- p1 + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
# p2 <- p2 + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))


gg <- (p1 / p2) +
    plot_layout(guides = "collect", widths = c(1, 1), heights = c(1.5,1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

ggsave(args[[2]], gg, width=22, height=22, units="cm")