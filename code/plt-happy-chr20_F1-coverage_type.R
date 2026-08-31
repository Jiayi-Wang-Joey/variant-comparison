#args <- list(list.files("~/omni/jiayiwang/data/variant-comparison/results/happy/",
#pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
    library(RColorBrewer)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt[, `:=`(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    METRIC.F1_Score        = 2 * TP / (2 * TP + FP + FN),
    TRUTH.TOTAL            = TP + FN
)]
.f <- \() {
    dt <- dt[!is.na(METRIC.F1_Score)]
    dt <- dt[!(grepl("longcallR", tool) & Type == "INDEL")]
    dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
    dt[,method:=paste(bamtype, aligner, tool, sep=">")]
    
    dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
    dt <- dt[platform == "MasSeq"]
    dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                            paste0("<span style='color:#54278f;'>", platform, "</span>"),
                            paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
    dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
    dt[cell_line=="HG002a", cell_line:="HG002"]
    dt[,Type:=factor(Type, levels = c("SNP", "INDEL"))]
    dt <- dt[aligner=="minimap2"]
    saveRDS(dt, "data/results/happy_chr20.rds")
    message(nrow(dt), " rows after filtering")
    cols <- c(
        "Clair3-RNA"   = "#A6CEE3",
        "DeepVariant"  = "#52AF43",
        "GATK"         = "#F06C45",
        "longcallR"    = "#B294C7",
        "longcallR-nn" = "#B15928",
        "isoLASER" = "#FDBF6F"
    )
    .p <- \(dt) {
        lb <- dt[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]),
                 by = .(coverage, Type, cell_line)]

        ggplot(dt, aes(coverage,
                       METRIC.F1_Score,
                       color = tool,
                       group = method)) +
            geom_point(alpha=0.8, size = 1.5) +
            geom_line(alpha=0.8, linewidth=0.8) +
            geom_text(data = lb,
                aes(x = coverage, y = 1.1, label = paste0(TRUTH.TOTAL)),
                inherit.aes = FALSE,
                angle = 30, hjust = 0.5, vjust = 1.3,
                size = 2, color = "grey30") +
            theme_classic() +
            facet_grid2(Type ~ cell_line, scales = "free") + 
            labs(
                x = "Coverage Cutoff (DP >= n)",
                y = "F1 Score",
                color = "Variant Caller",
            ) +
            scale_color_manual(values = cols) +
            theme(
                panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(
                    color = "black",
                    fill = NA,
                    linewidth = 0.8
                ),
                strip.background = element_rect(
                    fill = "white",
                    color = "black",
                    linewidth = 0.8),
                strip.text =  element_markdown(size=11),
                axis.line = element_line(color = "black", linewidth = 0.3),
                panel.spacing = unit(0, "lines"),
                panel.spacing.x = unit(0, "lines"),
                panel.spacing.y = unit(0, "lines"),
                axis.text.y = element_text(size = 7),
                axis.title.x = element_text(size = 11),
                axis.title.y = element_text(size = 11),
                legend.title = element_text(size = 11),
                aspect.ratio = 1
            )  
    }
    
    return(.p(dt))
    
}

gg <- .f()


ggsave(args[[2]], gg, width=27, height=13, units="cm")

#ggsave("omni/data/jiayiwang/variant-comparison/plts/happy-HG002.pdf", gg, width=22, height=12, units="cm")
