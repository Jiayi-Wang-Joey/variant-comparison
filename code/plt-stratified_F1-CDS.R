suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(patchwork)
    library(RColorBrewer)
})

cols <- c("Filter", "METRIC.F1_Score", "TRUTH.TOTAL", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR",tool) & Type == "INDEL")]
dt <- dt[Subtype == "*"]
dt <- dt[Subset=="refseq_cds"]
dt[,method:=paste(aligner, bamtype, tool, sep = ".")]



.f <- \(filter, title) {
    dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
    dt[,method:=paste(bamtype, aligner, tool, sep=">")]
    dt[, sample_lab := ifelse(grepl("MasSeq|IsoSeq", sample),
                              paste0("<span style='color:#8B4513;'>", sample, "</span>"),
                              paste0("<span style='color:#000000;'>", sample, "</span>"))]
    
    dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
    dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                            paste0("<span style='color:#54278f;'>", platform, "</span>"),
                            paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
    dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
    dt[cell_line == "HG002a", cell_line := "HG002"]
    dt <- dt[aligner=="minimap2"]
    snp <- dt[Type=="SNP"]
    idl <- dt[Type=="INDEL" & !grepl("longcallR", tool)]
    nk <- length(unique(dt$tool))
    cols <- c(
        "Clair3-RNA"   = "#A6CEE3",
        "DeepVariant"  = "#52AF43",
        "GATK"         = "#F06C45",
        "longcallR"    = "#B294C7",
        "longcallR-nn" = "#B15928",
        "isoLASER" = "#FDBF6F"
    )
    .p <- \(dt, title) {
        truth_labels <- dt[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(coverage, cell_line, platform)]
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
    p1 <- .p(snp, "SNP")
    p2 <- .p(idl, "INDEL") + theme(legend.position = "none")
    
    
    p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold")) 
}

gg <- .f("PASS", "F1 score - PASS and CDS only")


ggsave(args[[2]], gg, width=30, height=25, units="cm")