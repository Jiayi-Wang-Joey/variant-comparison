args <- list(list.files("/Volumes/jiayiwang/variant-comparison/results/SnpEff/",
pattern = "\\.vcf.gz", full.names = TRUE), "/Volumes/jiayiwang/variant-comparison/plts/snpeff_upset-alignment.pdf")

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(data.table)
    library(ComplexUpset)
    library(ggplot2)
    library(RColorBrewer)
    library(patchwork)
})
setwd("/Volumes/jiayiwang/variant-comparison/")
mmp <- readVcf("results/SnpEff/Clair3-RNA_origin_minimap2_HG005-MasSeq_5.vcf.gz")
pbm <- readVcf("results/SnpEff/Clair3-RNA_origin_pbmm2_HG005-MasSeq_5.vcf.gz")
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
    "Other"
) 
nk <- length(main_effects)
cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(nk),
                 unique(main_effects))
.p <- \(mmp, pbm, title, variant_type = "SNP") {
    
    .mkdt <- \(vcf, aligner, variant=variant_type) {
        #vcf <- vcf[geno(vcf)$BVT[,2]==variant,]
        dt <- as.data.table(geno(vcf)$BD)
        dt[,variant_id:=rownames(vcf)]
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
            #grepl("inframe_insertion|inframe_deletion", effects), "Inframe_indel",
            grepl("intron_variant", effects), "Intron",
            grepl("upstream_gene_variant|downstream_gene_variant", effects), "Upstream/Downstream",
            grepl("intergenic_region", effects), "Intergenic",
            grepl("non_coding_transcript", effects), "Noncoding Transcript",
            grepl("fusion", effects), "Fusion",
            default = "Other"
        )]
        
        dt <- na.omit(dt)
        dt[,type:=paste(aligner,res, sep="_")]
        dt
    }
    
    dt1 <- .mkdt(mmp, "minimap2")
    dt2 <- .mkdt(pbm, "pbmm2")
    dt <- rbind(dt1, dt2)
    
    types <- sort(unique(dt$type))
    
    m <- unique(dt[, .(variant_id, type, effect_main)])
    m_id <- unique(m[, .(variant_id, effect_main)])
    
    mem <- dcast(m, variant_id ~ type, fun.aggregate = length, value.var = "type")
    mem[, (types) := lapply(.SD, `>`, 0L), .SDcols = types]
    
    up <- merge(mem, m_id, by = "variant_id", all.x = TRUE)
    
    
    up2 <- up[!(minimap2_TP & pbmm2_TP)]
    ComplexUpset::upset(
        up2,
        types,
        set_sizes = FALSE,
        base_annotations = list(
            'Intersection size'=intersection_size(
                counts=FALSE,
                mapping=aes(fill=effect_main)
            ) + scale_fill_manual(values = cols) +
                theme(
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_line(),
                    panel.grid.minor.y = element_blank(),
                    axis.title.x = element_blank(),
                    legend.position = "bottom",
                    legend.box = "horizontal"
                ) +
                labs(fill="Variant Effect", x=NULL) +
                ggtitle(title)
                
        )
    )
}
caller <- c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "GATK")
sample <- c("HG004-MasSeq")
#sample <- c("HG004-Baylor-IsoSeq","HG005-Baylor-IsoSeq", "HG002-Baylor-IsoSeq",
#"HG004-MasSeq", "HG005-MasSeq","HG002-MasSeq")

params <- expand.grid(caller, sample)
ps <- lapply(seq_len(nrow(params)), \(i) {
    c <- params[i,1]
    s <- params[i,2]
    pbm <- readVcf(paste0("results/SnpEff/", c,"_origin_pbmm2_", s, "_5.vcf.gz"))
    mmp <- readVcf(paste0("results/SnpEff/", c,"_origin_minimap2_", s, "_5.vcf.gz"))
    #tt <- gsub("-Baylor", "", s)
    .p(mmp, pbm, title=c, variant_type = "SNP")
}) 
gg <- ps |> wrap_plots(ncol=5) + plot_layout(guide="collect")

ggsave("plts/upset_hg004-masseq.pdf", gg, width=55, height=12, units="cm")
#ggsave(args[[2]], ps, width=35, height=15, units="cm")
