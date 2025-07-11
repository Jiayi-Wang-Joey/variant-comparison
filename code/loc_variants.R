suppressPackageStartupMessages({
    library(txdbmaker)
    library(VariantAnnotation)
})

txdb <- makeTxDbFromGFF(args$gtf)
vcf <- readVcf(args$vcf)
chrs <- intersect(seqlevels(vcf), seqlevels(txdb))
vcf <- keepSeqlevels(vcf, chrs, pruning.mode = "coarse")
txdb <- keepSeqlevels(txdb, chrs, pruning.mode = "coarse")
loc <- locateVariants(vcf, txdb, AllVariants())
idx <- match(rownames(vcf),names(loc))
info(vcf)$LOCATION <- mcols(loc)$LOCATION[idx]
mcols(rowRanges(vcf))$LOCATION <- mcols(loc)$LOCATION[idx]

saveVcf(vcf, args$res)