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
# td <- td[
#     (tech %like% "IsoSeq" & status == "preprocessed") |  
#         (!tech %like% "IsoSeq" & status == "raw")            
# ]
td <- td[status=="raw"]
td[,error_rate:=error_rate*100]
td[, mapping_rate:=reads_mapped/raw_total_sequences]
td <- td[aligner=="minimap2"]
td <- td[,.(sample, status,  aligner, tech, reads_mapped,
            raw_total_sequences, mapping_rate, error_rate, 
            average_length)]
td[,raw_total_sequences:=raw_total_sequences/10e6]
td[,reads_mapped:=reads_mapped/10e6]
dd <- melt(
    td,
    id.vars = c("sample", "status", "aligner", "tech"),
    measure.vars = c("raw_total_sequences", "mapping_rate", "reads_mapped",
                     "error_rate", "average_length"),
    variable.name = "metric",
    value.name = "value"
)
saveRDS(dd, "data/results/qc.rds")
dd[sample=="HG002a", sample:="HG002"]
dd[tech=="Baylor-IsoSeq", tech:="IsoSeq"]
dd[, facet_group := interaction(metric, sample, drop = TRUE)]

metrics <- c(
    "raw_total_sequences",
    "reads_mapped",
    "mapping_rate",
    "average_length",
    "error_rate"
)
dd[, metric := factor(metric,
                       levels = metrics,
                       labels = c("Total reads (M)", "Mapped reads (M)", "Mapping rate (%)",
                                  "Average Read length", "Error rate (%)")
)]
dd <- dd[sample=="HG002"]
pal <- c(
    "dRNA002" = "#fff7bc",
    "cDNAxR09" = "#fee391",
    "dRNA004" = "#fec44f",
    "cDNAxR10" = "#d95f0e",
    "IsoSeq" = "#bcbddc",
    "MasSeq" = "#54278f"
)

dd[, col := pal[tech]]

aes <- list(
    geom_bar(stat = "identity", position = position_dodge()),
    scale_fill_identity(guide = "legend",
                        breaks = pal,
                        labels = names(pal)),
    facet_grid2(sample ~ metric, scales = "free", 
                axes="all", independent="all") ,
    theme_minimal(),
    theme(
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank()
    )
)
gg <- ggplot(dd, aes(reorder_within(tech, value, facet_group), 
                     value, fill=col)) + 
    aes +
    labs(x = NULL, y = "Value", fill = "Technology")


#nMetric <- length(unique(td$metric))
nSample <- length(unique(dd$sample))
ggsave(args[[2]], gg, width=24, height=4*nSample+0.5, units="cm")
