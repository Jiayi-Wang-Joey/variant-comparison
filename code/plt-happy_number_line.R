suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[
    Filter == "PASS" &
    !is.na(METRIC.F1_Score) &
    !(tool == "longcallR" & Type == "INDEL") &
    grepl("MasSeq|IsoSeq", sample)
]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]

td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "sample", "coverage"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]
p1 <- ggplot(snp, aes(coverage, Count, color=Category)) +
    geom_line(aes(linetype = aligner), alpha=0.6) +
    geom_point(aes(shape = aligner), alpha=0.6) +
    facet_grid2(tool ~ sample, scales="free", independent = "y") +
    theme_classic() +
    labs(y = "Total", x = "Variant caller") +
    scale_color_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    ggtitle("SNP")

p2 <- ggplot(idl, aes(coverage, Count, color=Category)) +
    geom_line(aes(linetype = aligner), alpha=0.6) +
    geom_point(aes(shape = aligner),  alpha=0.6) +
    facet_grid2(tool ~ sample, scales="free", independent = "y") +
    theme_classic() +
    labs(y = "Total", x = "Variant caller") +
    scale_color_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    ggtitle("INDEL")

gg <-   p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold")) 

ggsave(args[[2]], gg, width=30, height=28, units="cm")