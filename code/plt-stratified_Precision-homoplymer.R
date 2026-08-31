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
cols <- c("Filter", "METRIC.Precision", "TRUTH.TOTAL", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type", "Subset.Size")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS" & !is.na(METRIC.Precision)]
dt <- dt[!(grepl("longcallR",tool) & Type == "INDEL")]
dt <- dt[coverage==5 & Subtype=="*"]
dt <- dt[grepl("homopolymer", Subset)]
dt <- dt[aligner=="minimap2"]

dt[,method:=paste(aligner, bamtype, tool, sep = ".")]
dt$Subset <- recode(
    dt$Subset,
    "SimpleRepeat_homopolymer_4to6_AT_slop5"   = "4–6 AT",
    "SimpleRepeat_homopolymer_4to6_GC_slop5"   = "4–6 GC",
    "SimpleRepeat_homopolymer_4to6_slop5"      = "4–6",

    "SimpleRepeat_homopolymer_7to11_AT_slop5"  = "7–11 AT",
    "SimpleRepeat_homopolymer_7to11_GC_slop5"  = "7–11 GC",
    "SimpleRepeat_homopolymer_7to11_slop5"     = "7–11",

    "SimpleRepeat_homopolymer_ge12_AT_slop5"   = "≥12 AT",
    "SimpleRepeat_homopolymer_ge12_GC_slop5"   = "≥12 GC",
    "SimpleRepeat_homopolymer_ge12_slop5"      = "≥12",

    "SimpleRepeat_homopolymer_ge21_AT_slop5"   = "≥21 AT",
    "SimpleRepeat_homopolymer_ge21_GC_slop5"   = "≥21 GC",
    "SimpleRepeat_homopolymer_ge21_slop5"      = "≥21"
)

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt[cell_line == "HG002a", cell_line := "HG002"]
dt$Subset <- factor(
    dt$Subset,
    levels = c(
        "4–6", "4–6 AT", "4–6 GC",
        "7–11", "7–11 AT", "7–11 GC",
        "≥12", "≥12 AT", "≥12 GC",
        "≥21", "≥21 AT", "≥21 GC"
    )
)
dt_AT   <- dt[grepl("AT", Subset)]
dt_GC   <- dt[grepl("GC", Subset)]
dt_both <- dt[!grepl("AT|GC", Subset)]
cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER"     = "#FDBF6F"
)
.p <- \(dt) {
    dt <- dt[!(Type == "INDEL" & grepl("dRNA|cDNA", platform) & tool == "isoLASER")]
    snp <- dt[Type=="SNP"]
    idl <- dt[Type=="INDEL"]
    lb1 <- snp[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(Subset, cell_line, platform)]
    lb2 <- idl[, .(TRUTH.TOTAL = TRUTH.TOTAL[1]), by = .(Subset, cell_line, platform)]

    aes <- list(geom_point(alpha=0.6),
                geom_line(alpha = 0.6),
                facet_grid(cell_line ~ platform),
                theme_minimal(),
                labs(x = "Homopolymer length", y = "Precision", color = "Variant Caller"),
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
                      aspect.ratio = 1,
                     axis.text.x = element_text(angle = 45, size = 7,
                                                hjust = 1, vjust = 1)))

    p1 <- ggplot(snp, aes(Subset, METRIC.Precision,
                          col=tool,
                          group = method)) + aes +
        geom_text(data = lb1,
            aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
            inherit.aes = FALSE,
            angle = 30, hjust = 0.5, vjust = 1.3,
            size = 2, color = "grey30") +
        ggtitle("SNV")


    p2 <- ggplot(idl, aes(Subset, METRIC.Precision,
                          col=tool,
                          group = method)) + aes +
        geom_text(data = lb2,
            aes(x = Subset, y = 1.1, label = TRUTH.TOTAL),
            inherit.aes = FALSE,
            angle = 30, hjust = 0.5, vjust = 1.3,
            size = 2, color = "grey30") +
        ggtitle("INDEL") +
        theme(legend.position = "none")

    gg <- p1 + p2 + plot_layout(ncol = 1, guides = "collect") +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold"))
}

p1 <- .p(dt_both)
p2 <- .p(dt_AT)
p3 <- .p(dt_GC)

pdf(args[[2]], width = 9.4, height = 10.7)

p1
p2
p3

dev.off()
