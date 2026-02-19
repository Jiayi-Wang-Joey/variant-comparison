suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(patchwork)
    library(scales)
    library(dplyr)
    library(ggtext)
    library(RColorBrewer)
    library(tidytext)
})


#res <- lapply(args[[1]], fread, header=TRUE)
#res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
cols <- c("Filter", "METRIC.F1_Score", "coverage", "Subset", "Subtype",
          "sample", "aligner", "tool", "bamtype", "Type", "Subset.Size", "QUERY.TOTAL")
res <- lapply(args[[1]], function(f) {
    dt <- fread(f, header = TRUE)
    if (nrow(dt) > 0) dt[, ..cols] else NULL
})
res <- Filter(Negate(is.null), res)
dt <- rbindlist(res, use.names = TRUE)
dt <- dt[Filter=="PASS" & !is.na(METRIC.F1_Score)]
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
dt$Subset <- factor(
    dt$Subset,
    levels = c(
        "4–6", "4–6 AT", "4–6 GC",
        "7–11", "7–11 AT", "7–11 GC",
        "≥12", "≥12 AT", "≥12 GC",
        "≥21", "≥21 AT", "≥21 GC"
    )
)
dt[, text_col := fifelse(METRIC.F1_Score >= 0.8, "white", "black")]
dt_AT   <- dt[grepl("AT", Subset)]
dt_GC   <- dt[grepl("GC", Subset)]
dt_both <- dt[!grepl("AT|GC", Subset)]

.p <- \(dt) {
    snp <- dt[Type=="SNP"]
    idl <- dt[Type=="INDEL"]
    make_heat <- function(df) {
        df %>%
            group_by(cell_line, platform, Subset, tool, text_col) %>%
            summarise(
                METRIC.F1_Score = mean(METRIC.F1_Score, na.rm = TRUE),
                QUERY.TOTAL = sum(QUERY.TOTAL, na.rm = TRUE),   
                .groups = "drop"
            )
    }
    
    snp_h <- make_heat(snp)
    idl_h <- make_heat(idl)
    
    layers <- list(
        geom_tile(col = "white", linewidth = 0.1),
        geom_text(aes(label = QUERY.TOTAL,
                      color = text_col), 
                      size = 2.5),
        scale_fill_gradientn(
            "F1 score",
            colors = c("ivory", "gold", "red", "navy"),
            na.value = "lightgrey",
            limits = c(0, 1),
            n.breaks = 2
        ),
        scale_color_identity(),
        theme_minimal(),
        facet_grid2(cell_line ~ platform, scales = "free_y"),
        scale_y_reordered(sep = "___"),
        theme(
            plot.margin = margin(),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA),
            strip.text = element_markdown(),
            plot.tag = element_text(size = 9, face = "bold"),
            axis.text.x = element_text(angle = 45, size = 7,
                                       hjust = 1, vjust = 1)
        ),
        labs(x="Homopolymer length", y="Variant Caller"),
        guides(color = "none")
    )
    SEP <- "|" 
    p1 <- ggplot(
        snp_h,
        aes(
            Subset,
            reorder_within(tool, METRIC.F1_Score, cell_line, sep = SEP, desc = TRUE),
            fill = METRIC.F1_Score
        )
    ) +
        layers +
        scale_y_discrete(labels = function(x) sub("\\|.*$", "", x))
    
    p2 <- ggplot(
        idl_h,
        aes(
            Subset,
            reorder_within(tool, METRIC.F1_Score, cell_line, sep = SEP, desc = TRUE),
            fill = METRIC.F1_Score
        )
    ) +
        layers +
        scale_y_discrete(labels = function(x) sub("\\|.*$", "", x))
    
    gg <- p1 + p2 +
        plot_layout(ncol = 1, guides = "collect") +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold"))
}


p1 <- .p(dt_both)
p2 <- .p(dt_AT)
p3 <- .p(dt_GC)

pdf(args[[2]], width = 14, height = 13)

p1
p2
p3

dev.off()



write.table(dt, "data/results/homopolymers.csv")
    