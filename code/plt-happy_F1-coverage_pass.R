suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
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
    sample_levels <- dt[, .(platform = unique(platform)), by = sample][
        order(platform)
    ]$sample
    dt[, sample := factor(sample, levels = sample_levels)]
    dt[, sample := factor(sample, levels = sample_levels)]
    pb <- dt[grepl("IsoSeq|MasSeq", sample)]
    ont <- dt[!grepl("IsoSeq|MasSeq", sample)]
    p1 <- ggplot(ont, aes(coverage, 
                         METRIC.F1_Score, 
                         color = tool,
                         linetype = aligner,
                         shape = bamtype,
                         group = method)) +
        geom_point(alpha=0.6) +
        geom_line(alpha=0.6) + 
        facet_grid2(Type ~ sample, scales = "free") +
        theme_minimal() +
        labs(
            x = "Coverage Cutoff (DP >= n)",
            y = "F1 Score",
            color = "Variant Caller",
            linetype = "Aligner",
            shape = "Bam type"
        ) +
        scale_color_brewer(palette = "Set2") +
        theme(
            panel.border=element_rect(fill=NA)) +
        ggtitle("ONT")
    
    p2 <- ggplot(pb, aes(coverage, 
                                    METRIC.F1_Score, 
                                    color = tool,
                                    linetype = aligner,
                                    shape = bamtype,
                                    group = method)) +
        geom_point(alpha=0.6) +
        geom_line(alpha=0.6) + 
        facet_grid2(Type ~ sample, scales = "free") +
        theme_minimal() +
        labs(
            x = "Coverage Cutoff (DP >= n)",
            y = "F1 Score",
            color = "Variant Caller",
            linetype = "Aligner",
            shape = "Bam type"
        ) +
        scale_color_brewer(palette = "Set2") +
        theme(
            panel.border=element_rect(fill=NA)) +
        ggtitle("PacBio")
    
    p1 + p2 + plot_layout(ncol=1, guides = "collect") +
        plot_annotation(title=title, tag_levels = "a") &
        theme(
            plot.tag = element_text(face = "bold")  # Make tag bold
        )
}

gg <- .f("PASS", "F1 score - PASS only")

nSample <- length(unique(dt$sample))
ggsave(args[[2]], gg, width=7*nSample/2, height=25, units="cm")