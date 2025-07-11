# args <- list(c(
# list.files("/Volumes/jiayiwang/variant-comparison/results/align/preprocessed/quality/", full.names = TRUE),
# list.files("/Volumes/jiayiwang/variant-comparison/results/align/raw/quality/", full.names = TRUE)
# ))
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
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
td <- td[
    (tech %like% "IsoSeq" & status == "preprocessed") |  
        (!tech %like% "IsoSeq" & status == "raw")            
]
td[,error_rate:=error_rate*100]
td[, mapping_rate:=reads_mapped/raw_total_sequences]
td <- td[aligner=="minimap2"]
td <- td[,.(sample, status,  aligner, tech, 
            raw_total_sequences, mapping_rate, error_rate, 
            average_quality, average_length)]
td[,raw_total_sequences:=log10(raw_total_sequences)]
dd <- melt(
    td,
    id.vars = c("sample", "status", "aligner", "tech"),
    measure.vars = c("raw_total_sequences", "mapping_rate", 
                     "error_rate", "average_quality", "average_length"),
    variable.name = "metric",
    value.name = "value"
)


aes <- list(
    geom_bar(stat = "identity", position = position_dodge()),
    scale_fill_brewer(palette = "Paired"),
    facet_grid2(sample ~ metric, scales = "free", 
                axes="all", independent="all") ,
    theme_minimal(),
    theme(
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank()
    )
)


gg <- ggplot(dd, aes(reorder_within(tech, value, metric), 
                     value, fill=tech)) + 
    aes +
    labs(x = NULL, y = "Value", fill = "Technology")


#nMetric <- length(unique(td$metric))
nSample <- length(unique(td$sample))
ggsave(args[[2]], gg, width=25, height=4*nSample, units="cm")
