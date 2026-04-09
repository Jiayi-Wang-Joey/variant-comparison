# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")
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
dt <- dt[coverage==5]
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR", tool) & Type == "INDEL")]
dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[, sample_lab := ifelse(grepl("MasSeq|IsoSeq", sample),
                          paste0("<span style='color:#54278f;'>", sample, "</span>"),
                          paste0("<span style='color:#d95f0e;'>", sample, "</span>"))]

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[aligner=="minimap2"]
td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "cell_line", "platform"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]
p1 <- ggplot(snp, aes(reorder_within(tool,Count,platform), Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1
    ) +
    facet_grid2( cell_line ~ platform, scales="free") +
    scale_x_reordered() +
    theme_classic() +
    labs(y = "Total", x = "Variant caller") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    ggtitle("SNP")

p2 <- ggplot(idl, aes(tool, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1
    ) +
    facet_grid2( cell_line ~ platform, scales="free") + 
    theme_classic() +
    labs(y = "Total", x = "Variant caller") +
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

ggsave(args[[2]], gg, width=20, height=25, units="cm")

    
