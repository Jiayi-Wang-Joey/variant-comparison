# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(scales)
    library(RColorBrewer)
    library(ggtext)
    library(patchwork)
})

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(tool=="longcallR" & Type == "INDEL")]
dt[,coverage:=factor(coverage)]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[,aligner_tool:=paste(aligner, tool, sep = ".")]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
sample_levels <- dt[, .(platform = unique(platform)), by = sample][
    order(platform)
]$sample
dt[, sample := factor(sample, levels = sample_levels)]
dt <- dt[coverage==5]
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)
dt[, Type := factor(Type, c("SNP", "INDEL"))]
dt <- dt[aligner=="minimap2"]
dt[grepl("longcallR", tool) & Type=="INDEL", METRIC.Recall:=NA]
dt[grepl("longcallR", tool) & Type=="INDEL", METRIC.Precision:=NA]
# snp <- dt[Type=="SNP"]
# idl <- dt[Type=="INDEL" & !grepl("longcallR", tool)]

.p <- \(dt, title) {
    ggplot(dt, 
           aes(x = METRIC.Recall, 
               y = METRIC.Precision,
               group = method,
               col = tool,
               shape = cell_line)) +  
        geom_point(size=3, alpha = 0.8) +  
        facet_grid(Type ~ platform) +
        scale_color_manual(values = cols) +  
        labs(
            x = "Recall",
            y = "Precision",
            color = "Variant Caller",
            shape = "Cell line"
        ) +
        theme_classic() +
        scale_x_continuous(
            labels = scales::label_number(accuracy = 0.01)) +
        scale_y_continuous(
            labels = label_number(accuracy = 0.01)) +
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
            legend.position = "bottom"
        ) + coord_equal()
        
}

gg <- .p(dt) 

                
ggsave(args[[2]], gg, width=25, height=12, units="cm")
