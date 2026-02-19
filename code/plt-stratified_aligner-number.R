suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(tidytext)
    library(patchwork)
    library(RColorBrewer)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[
    coverage == 5 &
        Filter == "PASS" &
        !is.na(METRIC.F1_Score) &
        !(tool == "longcallR" & Type == "INDEL") &
        grepl("MasSeq|IsoSeq", sample) &
        Subtype == "*"
]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
keep_subsets <- c(
    "lowmappabilityall",
    # "nonunique_l100_m2_e1",
    # "nonunique_l250_m0_e0",	
    "notinlowmappabilityall"
)

dt <- dt[Subset %in% keep_subsets]

dt$Subset <- recode(dt$Subset,
                    "lowmappabilityall"="LM",
                    "notinlowmappabilityall"="NLM"
)
dt <- dt[grepl("MasSeq|IsoSeq", sample)]
cols <- colorRampPalette(brewer.pal(12, "Paired"))(4)
td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "sample", "Subset"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
td <- td[tool=="Clair3-RNA"]
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]
p1 <- ggplot(snp, aes(aligner, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1.5
    ) +
    facet_grid2(Subset ~ sample, scales="free") +
    scale_x_reordered() +
    theme_classic() +
    labs(y = "Total", x = "Variant Caller") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    ggtitle("SNP")

p2 <- ggplot(idl, aes(aligner, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1.5
    ) +
    facet_grid2(Subset ~ sample, scales="free") + 
    theme_classic() +
    labs(y = "Total", x = "Variant Caller") +
    scale_fill_brewer(palette = "Paired") +
    scale_x_reordered() +
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

ggsave(args[[2]], gg, width=25, height=28, units="cm")


