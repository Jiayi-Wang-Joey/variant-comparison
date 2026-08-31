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
        bamtype == "origin" &
        !is.na(METRIC.F1_Score) &
        !(grepl("longcallR", tool) & Type == "INDEL") &
        grepl("MasSeq|IsoSeq", sample) &
        !grepl("^HG002-", sample) &
        aligner %in% c("minimap2", "pbmm2")
]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "sample"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP", "METRIC.F1_Score"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
td[, sample:=gsub("-Baylor", "", sample)]
td[grepl("HG002a", sample), sample:=sub("HG002a", "HG002", sample)]
td <- unique(td)
saveRDS(td, "data/results/aligner.rds")
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]
# p1 <- ggplot(snp, aes(aligner, Count, fill=Category)) +
#     geom_col(stat = "identity") +
#     geom_text(
#         aes(label = Count),
#         position = position_stack(vjust = 0.5),
#         size = 1.5
#     ) +
#     facet_grid2(tool ~ sample, scales="free") +
#     scale_x_reordered() +
#     theme_classic() +
#     labs(y = "Total", x = "Variant caller") +
#     scale_fill_brewer(palette = "Paired") +
#     theme(
#         panel.border = element_rect(colour = "black", 
#                                     fill = NA, 
#                                     linewidth = 0.5),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text = element_markdown()
#     ) +
#     ggtitle("SNP")
p1 <- ggplot() +
    geom_col(
        data = snp[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, fill = Category)
    ) +
    geom_text(
        data = snp[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, label = Count, group = Category),
        position = position_stack(vjust = 0.5),
        #position = position_fill(vjust = 0.5),
        size = 1
    ) +
    geom_text(
        data = merge(
            snp[Category == "F1_Score", .(tool, sample, aligner, F1 = Count)],
            snp[Category %in% c("TP","FN","FP"),
                .(Total = sum(Count)),
                by = .(tool, sample, aligner)],
            by = c("tool","sample","aligner")
        ),
        aes(aligner, Total + 5000,
            label = sprintf("%.2f", F1)),
        size = 1.5,
        fontface = "bold"
    ) +
    facet_grid2(tool ~ sample, scales = "free") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    labs(y = "Total", x = "Variant caller") +
    ggtitle("SNV")

p2 <- ggplot() +
    geom_col(
        data = idl[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, fill = Category)
    ) +
    geom_text(
        data = idl[Category %in% c("TP","FN","FP")],
        aes(aligner, Count, label = Count, group=Category),
        position = position_stack(vjust = 0.5),
        size = 1
    ) +
    geom_text(
        data = merge(
            idl[Category == "F1_Score", .(tool, sample, aligner, F1 = Count)],
            idl[Category %in% c("TP","FN","FP"),
                .(Total = sum(Count)),
                by = .(tool, sample, aligner)],
            by = c("tool","sample","aligner")
        ),
        aes(aligner, Total + 1000,
            label = sprintf("%.2f", F1)),
        size = 1.5,
        fontface = "bold"
    ) +
    facet_grid2(tool ~ sample, scales = "free") +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_markdown()
    ) +
    labs(y = "Total", x = "Variant caller") +
    ggtitle("INDEL")

gg <-   p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold")) 

ggsave(args[[2]], gg, width=22, height=35, units="cm")


