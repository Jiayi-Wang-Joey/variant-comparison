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
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)
dt <- dt[!(grepl("longcallR", tool) & Type=="INDEL"),]

dt[,bamtype:=factor(bamtype, levels = c("origin", "transformed"))]
dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]


dt[,method:=paste(bamtype,tool,sep=".")]
dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
dt[,bamtype:=factor(bamtype, levels = c("origin", "transformed"))]
dt[,Type:=factor(Type, levels=c("SNP", "INDEL"))]
gg <- ggplot(dt, aes(coverage, 
               METRIC.F1_Score, 
               color = tool,
               group = method)) +
    geom_point(alpha=0.8, size = 1.5) +
    geom_line(alpha=0.8, inewidth=0.8, aes(linetype=bamtype)) + 
    theme_classic() +
    facet_grid2(Type ~ platform, scales = "free") + 
    labs(
        x = "Coverage Cutoff (DP >= n)",
        y = "F1 Score",
        color = "Variant Caller",
        linetype = "BAM Type"
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
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        #strip.text = element_text(size = 11),
        legend.text  = element_text(size = 9),
        aspect.ratio = 1
    ) 

ggsave(args[[2]], gg, width=30, height=12, units="cm")
#write.table(dt, "data/results/transformed_results.csv")

