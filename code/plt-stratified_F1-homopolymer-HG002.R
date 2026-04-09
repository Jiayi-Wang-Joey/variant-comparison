suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(RColorBrewer)
})


#res <- lapply(args[[1]], fread, header=TRUE)
#res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
cols <- c("Filter", "METRIC.F1_Score", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type", "Subset.Size")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS" & !is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR",tool) & Type == "INDEL")]
dt <- dt[coverage==5 & Subtype=="*"]
dt <- dt[grepl("homopolymer", Subset) & !grepl("AT|GC", Subset)]
dt <- dt[aligner=="minimap2"]
dt[,method:=paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(
    dt$Subset,
    "SimpleRepeat_homopolymer_4to6_slop5"      = "4–6",
    "SimpleRepeat_homopolymer_7to11_slop5"     = "7–11",
    "SimpleRepeat_homopolymer_ge12_slop5"      = "≥12",
    "SimpleRepeat_homopolymer_ge21_slop5"      = "≥21"
)

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[cell_line=="HG002"]
dt$Subset <- factor(
    dt$Subset,
    levels = c(
        "4–6", 
        "7–11", 
        "≥12", 
        "≥21"
    )
)
dt$Type <- factor(dt$Type, levels=c("SNP", "INDEL"))
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)
.p <- \(dt) {
    
    aes <- list(geom_point(alpha=0.6),
                geom_line(alpha = 0.6),
                facet_grid(Type ~ platform),
                theme_minimal(),
                labs(x = "Homopolymer length", y = "F1 Score", color = "Variant Caller"),
                scale_color_manual(values = cols),
                theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
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
                      #aspect.ratio = 1,
                      axis.text.x = element_text(angle = 45, size = 7,
                                                 hjust = 1, vjust = 1)))
    
    ggplot(dt, aes(Subset, METRIC.F1_Score, 
                          col=tool,
                          group = method)) + aes 
    }
gg <- .p(dt)

ggsave(args[[2]], gg, width=20, height=7, units="cm")

