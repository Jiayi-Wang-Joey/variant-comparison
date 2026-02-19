
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
})


res <- lapply(args[[1]], \(x) {
    dt <- fread(x)
    dt <- dt[, .(sample, bamtype, caller, phaser, 
                 aligner, chromosome, blockwise_hamming_rate)]
    #dt[, .(mean_switch_rate = mean(all_switch_rate)), 
    #by = .(sample, caller, phaser)]
})

dt <- rbindlist(res)
dt[, sample_lab := ifelse(grepl("MasSeq|IsoSeq", sample),
                          paste0("<span style='color:#8B4513;'>", 
                                 sample, "</span>"),
                          paste0("<span style='color:#000000;'>", 
                                 sample, "</span>"))]

gg <- ggplot(dt, aes(x = phaser, y = blockwise_hamming_rate, fill = phaser)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    facet_grid2(caller ~ sample_lab, scales="free_y", axes = "y") +
    theme_minimal(base_size = 8) +
    labs(y = "Switch Rate per Chromosome", x = "Phaser") +
    theme(
        panel.border=element_rect(fill=NA),
        strip.background = element_rect(fill = NA),
        axis.text.x = element_blank(),                    
        axis.ticks.x = element_blank(),
        strip.text = element_markdown())
scale_fill_brewer(palette = "Set1")  


ggsave(args[[2]], gg, width=30, height=8, units="cm")