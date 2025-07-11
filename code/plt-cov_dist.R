# args <- list(c(
#     list.files("/Volumes/jiayiwang/variant-comparison/results/align/preprocessed/coverage/", full.names = TRUE),
#     list.files("/Volumes/jiayiwang/variant-comparison/results/align/raw/coverage/", full.names = TRUE)
# ))

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

res <- lapply(args[[1]], fread)
dt <- rbindlist(res)
dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
dt[, sample := sub("^([^-]+)-.*$", "\\1", sample)]
dt <- dt[
    (tech %like% "IsoSeq" & status == "preprocessed") |  
        (!tech %like% "IsoSeq" & status == "raw")            
]
dt <- dt[aligner=="minimap2"]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", tech), "PacBio", "ONT")]

gg <- ggplot(dt, aes(x = depth, y = count, fill=platform)) +
    geom_bar(stat="identity") +
    scale_y_log10() +
    facet_grid(sample ~ tech) + 
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "black",
                                    fill = NA, linewidth = 0.5)) +
    scale_fill_manual(values = c("PacBio" = "#A0CED9", "ONT" = "tomato3"))
nSample <- length(unique(dt$sample))
nTech <- length(unique(dt$tech))

ggsave(args[[2]], gg, width=6*nTech, height=5*nSample, units="cm")

        
