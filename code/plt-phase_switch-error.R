# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/switch_error/",
# pattern = "\\.tsv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/phase_switch-error.pdf")

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

res <- lapply(args[[1]], \(x) {
    dt <- fread(x)
    dt <- dt[, .(sample, bamtype, caller, phaser, 
           aligner, chromosome, all_switch_rate)]
    dt[, .(mean_switch_rate = mean(all_switch_rate)), 
       by = .(sample, bamtype, caller, phaser, aligner)]
})

dt <- rbindlist(res)

ggplot(dt, aes(x = caller, y = mean_switch_rate, fill = phaser)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    facet_wrap(~sample+bamtype+aligner) +
    labs(y = "Mean Switch Rate", x = "Caller") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set1")

ggplot(dt, aes(x = caller, y = mean_switch_rate, fill = phaser)) +
    geom_violin(position = position_dodge(width = 0.9), trim = FALSE, alpha=0.5) +
    geom_jitter(
        aes(shape = bamtype, color = aligner),
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
        size = 2, alpha = 0.8
    ) +
    facet_grid(~sample) +
    theme_minimal() +
    labs(y = "Mean Switch Rate", x = "Caller") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 2),
           shape = guide_legend(order = 3))