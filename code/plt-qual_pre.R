# args <- list(c(
# list.files("/Volumes/jiayiwang/variant-comparison/results/align/preprocessed/quality/", full.names = TRUE),
# list.files("/Volumes/jiayiwang/variant-comparison/results/align/raw/quality/", full.names = TRUE)
# ))

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
})

res <- lapply(args[[1]], fread)
dt <- rbindlist(res)

# clean dt
dt$annotation <- NULL
dt[, tech := sub("^[^-]+-(.*)$", "\\1", sample)]
dt[, sample := sub("^([^-]+)-.*$", "\\1", sample)]
dt[, metric := sub(":$", "", metric)]
dt[, metric:=gsub(" ", "_", metric)]
td <- data.table::dcast(
    dt,
    sample + status + aligner + tech ~ metric,
    value.var = "value"
)
td <- td[grepl("IsoSeq", tech)]
td[,error_rate:=error_rate]
td[, mapping_rate:=reads_mapped/raw_total_sequences]
td <- td[aligner=="minimap2"]
td <- td[,.(sample, status,  aligner, tech, 
            raw_total_sequences, mapping_rate, error_rate, 
            average_quality, average_length)]

#td[,raw_total_sequences:=log10(raw_total_sequences)]
dd <- melt(td,
                id.vars = c("sample", "status", "aligner", "tech"),
                variable.name = "metric",
                value.name = "value"
)


gg <- ggplot(dd, aes(x = metric, y = value, fill = status)) +
    geom_col(position = position_dodge(width = 0.8),) +
    geom_text(aes(label = ifelse(metric == "error_rate",
                                 paste0(round(value * 100, 2), "%"),
                          ifelse(metric %in% c("raw_total_sequences", 
                                               "average_length"),
                                        format(round(value), big.mark = ",", 
                                               scientific = FALSE),
                                        format(round(value, 3), nsmall = 3)))),
              position = position_dodge(width = 0.8),
              vjust = 0, size = 3) +
    scale_fill_manual(values = c("preprocessed" = "#A0CED3", 
                                 "raw" = "tomato")) +
    facet_grid2(sample ~ metric, scales = "free", independent = "all") +
    labs(title = "Comparison of Preprocessed vs Raw Metrics",
         x = NULL, y = "Value", fill = "Status") +
    theme_minimal() +
    theme(
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank())

nSample <- length(unique(dd$sample))
ggsave(args[[2]], gg, width=33, height=7*nSample, units="cm")




