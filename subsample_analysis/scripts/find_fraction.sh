#!/usr/bin/env bash
# Binary search for the samtools-view subsample fraction that brings a source
# BAM's callable_bases (at a given coverage cutoff) closest to a target value.
#
# Each iteration mirrors workflow/rules/variant_benchmark.smk:generate_bed
# exactly (mosdepth -> DP filter -> bedtools merge -> exon intersect ->
# confidence-region intersect), so the candidate callable_bases number is
# computed the same way the rest of the pipeline computes it.
set -euo pipefail

usage() {
    echo "usage: $0 <bam> <coverage> <gtf> <confidence_bed> <target_tsv> <mosdepth_img> <seed> <max_iter> <tol> <threads> <tmpdir> <out_fraction> <out_log>" >&2
    exit 1
}

[ "$#" -eq 13 ] || usage

BAM="$1"
COVERAGE="$2"
GTF="$3"
CONFIDENCE="$4"
TARGET_TSV="$5"
MOSDEPTH_IMG="$6"
SEED="$7"
MAX_ITER="$8"
TOL="$9"
THREADS="${10}"
TMPDIR="${11}"
OUT_FRACTION="${12}"
OUT_LOG="${13}"

mkdir -p "$TMPDIR" "$(dirname "$OUT_FRACTION")" "$(dirname "$OUT_LOG")"

TARGET=$(awk -F'\t' '$4=="callable_bases"{print $5}' "$TARGET_TSV")
if [ -z "$TARGET" ]; then
    echo "ERROR: could not find callable_bases in $TARGET_TSV" >&2
    exit 1
fi

echo "iter	fraction	callable_bases	target	diff	rel_diff" > "$OUT_LOG"

EXON_BED="$TMPDIR/exons.bed"
grep -v "^#" "$GTF" | awk 'BEGIN{OFS="\t"} $3=="exon"' > "$EXON_BED"

low="0.0"
high="1.0"
best_frac="1.0"
best_diff=""

for i in $(seq 1 "$MAX_ITER"); do
    frac=$(awk -v l="$low" -v h="$high" 'BEGIN{printf "%.4f", (l+h)/2}')
    subsample_param=$(awk -v s="$SEED" -v f="$frac" 'BEGIN{printf "%.4f", s+f}')

    iter_bam="$TMPDIR/iter${i}.bam"
    prefix="$TMPDIR/iter${i}"

    samtools view -s "$subsample_param" -@ "$THREADS" -b "$BAM" -o "$iter_bam"
    samtools index "$iter_bam"

    singularity exec \
        --bind "$(pwd -P)":"$(pwd -P)" \
        --pwd "$(pwd -P)" \
        "$MOSDEPTH_IMG" mosdepth \
        --threads "$THREADS" \
        "$prefix" \
        "$iter_bam"

    callable=$(
        gzip -fdc "${prefix}.per-base.bed.gz" \
            | awk -v min_cov="$COVERAGE" '$4 >= min_cov' \
            | bedtools merge -d 1 -c 4 -o mean -i - \
            | bedtools intersect -a - -b "$EXON_BED" \
            | bedtools intersect -a - -b "$CONFIDENCE" \
            | awk 'BEGIN{s=0} {s += $3-$2} END{print s+0}'
    )

    diff=$(awk -v c="$callable" -v t="$TARGET" 'BEGIN{d=c-t; if (d<0) d=-d; print d}')
    rel_diff=$(awk -v d="$diff" -v t="$TARGET" 'BEGIN{printf "%.6f", d/t}')

    echo "${i}	${frac}	${callable}	${TARGET}	${diff}	${rel_diff}" >> "$OUT_LOG"

    if [ -z "$best_diff" ] || awk -v d="$diff" -v b="$best_diff" 'BEGIN{exit !(d<b)}'; then
        best_diff="$diff"
        best_frac="$frac"
    fi

    rm -f "$iter_bam" "${iter_bam}.bai" \
        "${prefix}.per-base.bed.gz" "${prefix}.per-base.bed.gz.csi" \
        "${prefix}.mosdepth.global.dist.txt" "${prefix}.mosdepth.summary.txt" \
        "${prefix}.mosdepth.region.dist.txt"

    if awk -v c="$callable" -v t="$TARGET" 'BEGIN{exit !(c>t)}'; then
        high="$frac"
    else
        low="$frac"
    fi

    if awk -v rd="$rel_diff" -v tol="$TOL" 'BEGIN{exit !(rd<=tol)}'; then
        break
    fi
done

echo "$best_frac" > "$OUT_FRACTION"
rm -f "$EXON_BED"

echo "Chosen fraction: $best_frac (best diff seen: $best_diff, target: $TARGET)" >&2
