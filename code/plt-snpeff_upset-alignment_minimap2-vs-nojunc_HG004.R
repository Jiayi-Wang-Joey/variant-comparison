suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
    library(ComplexUpset)
    library(ggplot2)
    library(RColorBrewer)
    library(patchwork)
    library(grid)
})
setwd("/data/jiayiwang/variant-comparison/")

main_effects <- c(
    "LoF",
    "Missense",
    "Synonymous",
    "Splice Region",
    "UTR",
    "Intron",
    "Upstream/Downstream",
    "Intergenic",
    "Noncoding Transcript",
    "Fusion",
    "Other"
)
nk <- length(main_effects)
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk), main_effects)

.p <- \(mmp, njc, title, variant_type = "SNP", legend = FALSE) {

    .mkdt <- \(vcf, aligner, variant = variant_type) {
        dt <- as.data.table(geno(vcf)$BD)
        dt[, variant_id := rownames(vcf)]
        dt[, res := fcase(
            TRUTH == "TP", "TP",
            QUERY == "FP" & TRUTH == ".", "FP",
            TRUTH == "FN" & QUERY == "FP", "FN",
            default = NA_character_
        )]
        dt$effects <- sub(
            "^[^|]*\\|([^|]+)\\|.*$",
            "\\1",
            vapply(info(vcf)$ANN, `[`, character(1), 1)
        )
        dt[, effect_main := fcase(
            grepl("frameshift_variant|stop_gained|stop_lost|start_lost|transcript_ablation", effects), "LoF",
            grepl("missense_variant", effects), "Missense",
            grepl("synonymous_variant|stop_retained_variant|start_retained_variant", effects), "Synonymous",
            grepl("splice_acceptor_variant|splice_donor_variant|splice_region_variant", effects), "Splice Region",
            grepl("3_prime_UTR|5_prime_UTR", effects), "UTR",
            grepl("intron_variant", effects), "Intron",
            grepl("upstream_gene_variant|downstream_gene_variant", effects), "Upstream/Downstream",
            grepl("intergenic_region", effects), "Intergenic",
            grepl("non_coding_transcript", effects), "Noncoding Transcript",
            grepl("fusion", effects), "Fusion",
            default = "Other"
        )]

        dt <- na.omit(dt)
        dt[, type := paste(aligner, res, sep = "_")]
        dt
    }

    dt1 <- .mkdt(mmp, "minimap2")
    dt2 <- .mkdt(njc, "nojunc")
    dt <- rbind(dt1, dt2)

    m <- unique(dt[, .(variant_id, type, effect_main)])
    m_id <- unique(m[, .(variant_id, effect_main)])

    types_all <- sort(unique(dt$type))
    mem <- dcast(m, variant_id ~ type, fun.aggregate = length, value.var = "type")
    mem[, (types_all) := lapply(.SD, `>`, 0L), .SDcols = types_all]

    up <- merge(mem, m_id, by = "variant_id", all.x = TRUE)

    up2 <- up[(minimap2_FP | nojunc_FP) &
                  !(minimap2_TP | minimap2_FN | nojunc_TP | nojunc_FN)]
    types <- c("minimap2_FP", "nojunc_FP")

    ComplexUpset::upset(
        up2,
        types,
        set_sizes = FALSE,
        base_annotations = list(
            "Intersection size" = intersection_size(
                counts = FALSE,
                mapping = aes(fill = effect_main)
            ) + scale_fill_manual(values = cols, drop = FALSE) +
                theme(
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_line(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = if (legend) "right" else "none"
                ) +
                labs(fill = "Variant Effect") +
                ggtitle(title)
        )
    )
}

caller <- c("DeepVariant", "Clair3-RNA", "longcallR", "longcallR-nn", "isoLASER")
sample <- c("MasSeq", "Baylor-IsoSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004")
excluded_sample <- "cDNAxR10"
included_samples <- setdiff(sample, excluded_sample)

# one row per dataset: a rotated dataset label on the left + one upset plot per caller
.mkrow <- \(s, legend = FALSE) {
    ps <- lapply(seq_along(caller), \(i) {
        c <- caller[i]
        njc <- readVcf(paste0("results/SnpEff/", c, "_origin_minimap2-nojunc_HG004-", s, "_5.vcf.gz"))
        mmp <- readVcf(paste0("results/SnpEff/", c, "_origin_minimap2_HG004-", s, "_5.vcf.gz"))
        .p(mmp, njc, title = c, variant_type = "SNP", legend = legend && i == 1)
    })
    label <- wrap_elements(grid::textGrob(
        paste0("HG004-", s), rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")
    ))
    cat("done:", s, "\n")
    (label | wrap_plots(ps, ncol = length(caller))) + plot_layout(widths = c(0.025, 1))
}

row_w <- 36 / 2.54
extra_w <- 6 / 2.54
page_w <- row_w + extra_w
row_h <- 10 / 2.54

out <- "plts/upset_FPs_minimap2-vs-nojunc_HG004.pdf"
tmp_main <- tempfile(fileext = ".pdf")
tmp_excl <- tempfile(fileext = ".pdf")

# page 1: all included datasets stacked in a single column, sharing one legend
pdf(tmp_main, width = page_w, height = row_h * length(included_samples))
rows <- lapply(included_samples, \(s) .mkrow(s, legend = s == included_samples[1]))
print(wrap_plots(rows, ncol = 1) + plot_layout(guides = "collect"))
dev.off()

# page 2: the excluded dataset on its own page, same legend styling
pdf(tmp_excl, width = page_w, height = row_h)
row_excl <- .mkrow(excluded_sample, legend = TRUE)
print(wrap_plots(list(row_excl), ncol = 1) + plot_layout(guides = "collect"))
dev.off()

status <- system2(
    "gs",
    c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", paste0("-o", out), tmp_main, tmp_excl)
)
if (status != 0) stop("failed to merge PDF pages with ghostscript")
file.remove(tmp_main, tmp_excl)
cat("done:", out, "\n")
