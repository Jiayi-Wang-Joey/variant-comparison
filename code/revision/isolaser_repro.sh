#!/usr/bin/env bash
# isoLASER at the cutoffs of Quinones-Valdez et al. 2025 (Fig. 2a), HG002 Iso-Seq
# Usage: isolaser_repro.sh <outdir>
set -euo pipefail
cd /data/jiayiwang/variant-comparison
O=${1:-results/revision/isolaser_repro}
mkdir -p $O
BIN=.snakemake/conda/f8d0742757163c52f72d7a18cbbde687_/bin
export PATH=$BIN:$PATH TMPDIR=/data/jiayiwang/tmp
mkdir -p $TMPDIR
S=HG002a-Baylor-IsoSeq
GTF=/home/jiayiwang/data/reference/gencode.v47.chr_patch_hapl_scaff.annotation.sorted.gtf
CONF=data/truth/beds/HG002.bed
TRUTH=data/truth/HG002/truth.vcf.gz
REF=data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta
CUTS="20 31 100 310"

grep -v "^#" $GTF | awk '$3=="exon"' > $O/exons.gtf
for c in $CUTS; do
    gzip -dc results/beds/raw/$S/minimap2_100.per-base.bed.gz | awk -v m=$c '$4>=m' \
        | bedtools merge -d 1 -c 4 -o mean -i - \
        | bedtools intersect -a - -b $O/exons.gtf | bedtools intersect -a - -b $CONF > $O/raw_$c.bed
done
rm -f $O/exons.gtf

Q=$O/isoLASER.query.vcf.gz
bcftools view -f PASS -Oz -o $Q results/variant/isoLASER_filtered/origin_minimap2_$S.vcf.gz
bcftools index -t -f $Q
for c in $CUTS; do
    HGREF=$REF singularity exec --bind $(pwd -P):$(pwd -P) --bind $TMPDIR:$TMPDIR --pwd $(pwd -P) \
        /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/hap.py $TRUTH $Q -r $REF -f $CONF -T $O/raw_$c.bed \
        -o $O/isoLASER_$c --engine=vcfeval --scratch-prefix $TMPDIR --threads 8 > $O/isoLASER_$c.log 2>&1
done
echo DONE > $O/done.txt
