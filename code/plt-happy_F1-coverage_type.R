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

.f <- \(filter, title) {
    dt <- dt[Filter==filter]
    dt <- dt[!is.na(METRIC.F1_Score)]
    dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
    dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
    dt[,method:=paste(bamtype, aligner, tool, sep=">")]
    
    dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
    dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                              paste0("<span style='color:#54278f;'>", platform, "</span>"),
                              paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
    dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
    dt <- dt[aligner=="minimap2"]
    dt[cell_line=="HG002a", cell_line:="HG002"]
    saveRDS(dt, "data/results/happy.csv")
    dt <- dt[cell_line != "HG002"]
    snp <- dt[Type=="SNP"]
    idl <- dt[Type=="INDEL" & !grepl("longcallR", tool)]
    idl <- idl[!(grepl("dRNA|cDNA", platform) & tool == "isoLASER")]
    cols <- c(
        "Clair3-RNA"   = "#A6CEE3",
        "DeepVariant"  = "#52AF43",
        "GATK"         = "#F06C45",
        "longcallR"    = "#B294C7",
        "longcallR-nn" = "#B15928",
        "isoLASER" = "#FDBF6F"
    )
    .p <- \(dt, title, truth_labels) {
        ggplot(dt, aes(coverage, 
                        METRIC.F1_Score, 
                        color = tool,
                        group = method)) +
            geom_point(alpha=0.8, size = 1.5) +
            geom_line(alpha=0.8, linewidth=0.8) + 
            geom_text(data = truth_labels,
                aes(x = coverage, y = 1.1, label = paste0(TRUTH.TOTAL)),
                inherit.aes = FALSE,
                angle = 30, hjust = 0.5, vjust = 1.3,
                size = 2, color = "grey30") +
            theme_classic() +
            facet_grid2( cell_line ~ platform, scales = "free") + 
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
            )  +
            ggtitle(title)
    }
    lb1 <- snp[Type=="SNP", .(TRUTH.TOTAL = TRUTH.TOTAL[1]),
                   by = .(coverage, cell_line, platform)]
    lb2 <- idl[Type=="INDEL", .(TRUTH.TOTAL = TRUTH.TOTAL[1]),
                by = .(coverage, cell_line, platform)]
    p1 <- .p(snp, "SNP", lb1) #+ theme(legend.position = "none")
    p2 <- .p(idl, "INDEL", lb2) + theme(legend.position = "none")
    
    
    p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold")) 
}

gg <- .f("PASS", "F1 score - PASS only")


ggsave(args[[2]], gg, width=30, height=25, units="cm")

#ggsave("omni/data/jiayiwang/variant-comparison/plts/happy-HG002.pdf", gg, width=22, height=12, units="cm")
