
# args <- list(list.files("~/omni/data/jiayiwang/variant-comparison/results/align/raw/supplementary_stats/",
# pattern = "\\.tsv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
    library(tidyverse)
})

res <- lapply(args[[1]], fread)
dt <- rbindlist(res)
saveRDS(dt, "data/supp_rate.rds")
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[status=="raw"]

td <- dt %>%
  pivot_wider(
    id_cols = c(metric, platform, cell_line, sample),
    names_from = aligner,
    values_from = value
  )

pal <- c(
    "dRNA002" = "#fff7bc",
    "cDNAxR09" = "#fee391",
    "dRNA004" = "#fec44f",
    "cDNAxR10" = "#d95f0e",
    "IsoSeq" = "#bcbddc",
    "MasSeq" = "#54278f"
)
setDT(td)
td <- td[metric=="supplementary_read_rate"]


gg <- ggplot(td, aes(minimap2, pbmm2, col = platform)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = pal) +
  labs(x = "minimap2", y = "pbmm2", color = "Platform") +
  theme_classic() +
  facet_wrap(~cell_line, scales = "free") 

ggsave(args[[2]], gg, width = 6, height = 4)