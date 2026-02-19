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
        coverage == 5 &
        Filter == "PASS" &
        !is.na(METRIC.F1_Score) &
        aligner == "minimap2" &
        !(grepl("longcallR",tool) & Type == "INDEL") ]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
td <- melt(
    dt,
    id.vars = c("Type", "aligner", "tool", "platform", "bamtype"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]
p1 <- ggplot(snp, aes(bamtype, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1.5
    ) +
    facet_grid2(tool ~ platform, scales="free") +
    scale_x_reordered() +
    theme_classic() +
    labs(y = "Total", x = "BAM type") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown(),
    ) +
    ggtitle("SNP")

p2 <- ggplot(idl, aes(bamtype, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1.5
    ) +
    facet_grid2(tool ~ platform, scales="free") +
    scale_x_reordered() +
    theme_classic() +
    labs(y = "Total", x = "BAM type") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown(),
    ) +
    ggtitle("SNP")

gg <-   p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold")) 

ggsave(args[[2]], gg, width=28, height=28, units="cm")
