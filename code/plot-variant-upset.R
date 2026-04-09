suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
    library(ComplexUpset)
    library(ggplot2)
    library(RColorBrewer)
    library(patchwork)
})
setwd("/Volumes/jiayiwang/variant-comparison/")
tech <- c("MasSeq", "cDNAxR10", "dRNA004")

get_dt <- function(vcf, caller) {
    dt <- as.data.table(geno(vcf)$BD)
    dt[, variant_id := rownames(vcf)]
    
    dt[, res := fcase(
        QUERY == "FP", "FP",
        TRUTH == "FN", "FN",
        TRUTH == "TP" | QUERY == "TP", "TP",
        TRUTH == "N", NA_character_,
        default = NA_character_
    )]
    
    dt[, caller := caller]
    dt[, .(variant_id, caller, res)]
}

plot_upset <- function(tech) {
    clair3 <- readVcf(paste0("results/happy/Clair3-RNA_origin_minimap2_HG004-", tech, "_5.vcf.gz"))
    dv <- readVcf(paste0("results/happy/DeepVariant_origin_minimap2_HG004-", tech, "_5.vcf.gz"))
    longcallR <- readVcf(paste0("results/happy/longcallR_origin_minimap2_HG004-", tech, "_5.vcf.gz"))
    longcallRnn <- readVcf(paste0("results/happy/longcallR-nn_origin_minimap2_HG004-", tech, "_5.vcf.gz"))
    
    long_dt <- rbindlist(list(
        get_dt(clair3, "Clair3-RNA"),
        get_dt(dv, "DeepVariant"),
        get_dt(longcallR, "longcallR"),
        get_dt(longcallRnn, "longcallR-nn")
    ), use.names = TRUE)
    
    long_dt <- unique(long_dt, by = c("variant_id", "caller"))
    
    presence_dt <- dcast(
        long_dt,
        variant_id ~ caller,
        value.var = "caller",
        fun.aggregate = length
    )
    
    callers <- c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn")
    
    for (cc in callers) {
        presence_dt[, (cc) := as.integer(get(cc) > 0)]
    }
    
    ## one label per variant for coloring
    ## priority: FP > FN > TP
    res_dt <- long_dt[!is.na(res), .(
        res = fcase(
            any(res == "FP"), "FP",
            any(res == "FN"), "FN",
            any(res == "TP"), "TP",
            default = "Other"
        )
    ), by = variant_id]
    
    plot_dt <- merge(presence_dt, res_dt, by = "variant_id", all.x = TRUE)
    plot_dt[is.na(res), res := "Other"]
    
    paired_cols <- brewer.pal(12, "Paired")
    fill_cols <- c(
        TP = paired_cols[2],
        FP = paired_cols[6],
        FN = paired_cols[4],
        Other = "grey70"
    )
    
    upset(
        plot_dt,
        intersect = callers,
        name = "Variants",
        base_annotations = list(
            "Intersection size" = intersection_size(
                aes(fill = res)
            )
        ),
        width_ratio = 0.18,
        min_size = 1
    ) +
        scale_fill_brewer(palette = "Paired") +
        ggtitle(paste("HG004 -", tech)) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}
