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

dt <- dt[aligner=="minimap2"]
dt[,method:=paste(aligner, bamtype, tool, sep = ".")]

keep_subsets <- c(
    "SimpleRepeat_diTR_10to49_slop5",
    "SimpleRepeat_diTR_50to149_slop5",
    "SimpleRepeat_diTR_ge150_slop5",
    "SimpleRepeat_triTR_14to49_slop5",	
    "SimpleRepeat_triTR_50to149_slop5",	
    "SimpleRepeat_triTR_ge150_slop5",
    "SimpleRepeat_quadTR_19to49_slop5", 
    "SimpleRepeat_quadTR_50to149_slop5",
    "SimpleRepeat_quadTR_ge150_slop5"
)
dt <- dt[Subset %in% keep_subsets]
dt$Subset <- recode(dt$Subset,
                     "SimpleRepeat_diTR_10to49_slop5"="diTR_10-49",
                     "SimpleRepeat_diTR_50to149_slop5"="diTR_50-149",
                     "SimpleRepeat_diTR_ge150_slop5"="diTR_>150",
                     "SimpleRepeat_triTR_14to49_slop5"="triTR_14-49",	
                     "SimpleRepeat_triTR_50to149_slop5"="triTR_50-149",	
                     "SimpleRepeat_triTR_ge150_slop5"="triTR_>150",
                     "SimpleRepeat_quadTR_19to49_slop5"="quadTR_19-49", 
                     "SimpleRepeat_quadTR_50to149_slop5"="quadTR_50-149",
                     "SimpleRepeat_quadTR_ge150_slop5"="quadTR_>150"
)

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
lvls <- c(
    "diTR_10-49", "diTR_50-149","diTR_>150",
    "triTR_14-49", "triTR_50-149","triTR_>150",
    "quadTR_19-49", "quadTR_50-149", "quadTR_>150"
)

dt$Subset <- factor(
    dt$Subset,
    levels = lvls
)
dt[, text_col := fifelse(METRIC.F1_Score >= 0.8, "white", "black")]

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


p1 <- .p(dt)


pdf(args[[2]], width = 14, height = 13)

p1


dev.off()