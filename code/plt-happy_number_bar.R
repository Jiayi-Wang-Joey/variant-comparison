# args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/happy/",
# pattern = "\\.summary.csv$", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/happy_F1-bar.pdf")
suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
    library(patchwork)
    library(ggtext)
})

# y-positions for stacked-bar count labels, nudged apart to avoid overlap
add_label_pos <- function(dt, by, scale_by = by, min_gap_frac = 0.075) {
    dt <- copy(dt)
    setorderv(dt, c(by, "Category"))
    dt[, ymax := rev(cumsum(rev(Count))), by = by]
    dt[, ymin := ymax - Count]
    dt[, ymid := (ymin + ymax) / 2]
    dt[, total := sum(Count), by = by]
    dt[, scale_ref := max(total), by = scale_by]
    dt[, label_y := {
        o <- order(ymid)
        y <- ymid[o]
        mg <- scale_ref[1] * min_gap_frac
        if (length(y) > 1) {
            for (i in 2:length(y)) if (y[i] - y[i - 1] < mg) y[i] <- y[i - 1] + mg
        }
        out <- numeric(length(y))
        out[o] <- y
        out
    }, by = by]
    dt
}

res <- lapply(args[[1]], fread, header=TRUE)
res <- res[!vapply(res, \(.) nrow(.)==0, logical(1))]
dt <- rbindlist(res)
dt <- dt[coverage==5]
dt <- dt[Filter=="PASS"]
dt <- dt[!is.na(METRIC.F1_Score)]
dt <- dt[!(grepl("longcallR|NanoTS", tool) & Type == "INDEL")]
dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
dt[,method:=paste(bamtype, aligner, tool, sep=">")]
dt[, sample_lab := ifelse(grepl("MasSeq|IsoSeq", sample),
                          paste0("<span style='color:#54278f;'>", sample, "</span>"),
                          paste0("<span style='color:#d95f0e;'>", sample, "</span>"))]

dt[, platform := factor(sapply(strsplit(sample, "-"), tail, 1))]
dt[, platform := ifelse(grepl("MasSeq|IsoSeq", sample),
                        paste0("<span style='color:#54278f;'>", platform, "</span>"),
                        paste0("<span style='color:#d95f0e;'>", platform, "</span>"))]
dt[, cell_line := factor(sapply(strsplit(sample, "-"), head, 1))]
dt <- dt[aligner=="minimap2"]
dt <- dt[cell_line=="HG002a", cell_line:="HG002"]
# isoLASER indels: PacBio only
dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & grepl("dRNA|cDNA", platform))]
td <- melt(
    dt,
    id.vars = c("Type", "Filter", "aligner", "tool", "cell_line", "platform"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
td[, Category := tstrsplit(Category, "\\.", keep = 2)]
td <- unique(td)
snp <- td[Type=="SNP"]
idl <- td[Type=="INDEL"]

snp_lab <- add_label_pos(snp, by = c("tool", "cell_line", "platform"), scale_by = "cell_line")
idl_lab <- add_label_pos(idl, by = c("tool", "cell_line", "platform"), scale_by = "cell_line")

p1 <- ggplot(snp, aes(reorder_within(tool,Count,platform), Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        data = snp_lab,
        aes(x = reorder_within(tool, Count, platform), y = label_y, label = Count),
        inherit.aes = FALSE, size = 2.2
    ) +
    facet_grid2( cell_line ~ platform, scales="free") +
    scale_x_reordered() +
    theme_classic(base_size = 11) +
    labs(y = "Total", x = "Variant caller") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_markdown(size = 8)
    ) +
    ggtitle("SNV")

p2 <- ggplot(idl, aes(tool, Count, fill=Category)) +
    geom_col(stat = "identity") +
    geom_text(
        data = idl_lab,
        aes(x = tool, y = label_y, label = Count),
        inherit.aes = FALSE, size = 2.2
    ) +
    facet_grid2( cell_line ~ platform, scales="free") +
    theme_classic(base_size = 11) +
    labs(y = "Total", x = "Variant caller") +
    scale_fill_brewer(palette = "Paired") +
    scale_x_reordered() +
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_markdown(size = 8)
    ) +
    ggtitle("INDEL")

write.table(dt, "data/results/number_calls.csv")
pdf(args[[2]], width = 26 / 2.54, height = 18 / 2.54)
print(p1)
print(p2)
dev.off()

    
