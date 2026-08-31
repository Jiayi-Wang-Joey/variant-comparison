#!/usr/bin/env python3
# summarize_phase_blocks.py
# Combine three `whatshap stats` outputs into a single summary row with
# phase-block N50:
#   --tsv         per-chromosome + ALL summary (incl. fraction of
#                 heterozygous variants phased)
#   --block-list  one row per phase block/singleton (sample, chromosome,
#                 phase_set, from, to, variants)
#   --gtf         each block as a "gene", split into multiple "exon" lines
#                 when blocks are geometrically interleaved/nested. This is
#                 needed because summing (to - from) directly from
#                 --block-list double-counts bp for overlapping blocks;
#                 summing exon lengths per gene_id gives whatshap's own,
#                 non-overlapping block length (verified against bp_per_
#                 block_sum in --tsv on real data).
#
# Only blocks with >=2 variants count towards N50 (singletons carry no
# linkage information), matching whatshap's own "Blocks" vs "Singletons"
# distinction in --tsv.
#
# Usage: python3 summarize_phase_blocks.py --blocks blocks.tsv --gtf blocks.gtf
#            --stats stats.tsv --sample S --caller C --phaser P --type T
#            --output out.tsv

import argparse
import csv
import re
from collections import defaultdict


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def read_gtf_block_lengths(path):
    lengths = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            chrom, _, _, start, end, *_, attrs = line.rstrip("\n").split("\t")
            gene_id = re.search(r'gene_id "([^"]+)"', attrs).group(1)
            lengths[(chrom, gene_id)] += int(end) - int(start)
    return lengths


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    if total == 0:
        return 0
    cum = 0
    for length in lengths:
        cum += length
        if cum >= total / 2:
            return length
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--phaser", required=True)
    ap.add_argument("--type", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    bp_by_block = read_gtf_block_lengths(args.gtf)

    blocks = [b for b in read_tsv(args.blocks) if int(b["variants"]) >= 2]
    lengths_bp = [bp_by_block[(b["chromosome"], b["phase_set"])] for b in blocks]
    lengths_variants = [int(b["variants"]) for b in blocks]

    stats = read_tsv(args.stats)
    all_row = next((r for r in stats if r["chromosome"] == "ALL"), {})

    row = {
        "sample": args.sample,
        "caller": args.caller,
        "phaser": args.phaser,
        "type": args.type,
        "n_variants": all_row.get("variants", "NA"),
        "n_phased": all_row.get("phased", "NA"),
        "n_unphased": all_row.get("unphased", "NA"),
        "n_singletons": all_row.get("singletons", "NA"),
        "n_blocks": len(blocks),
        "block_n50_bp": n50(lengths_bp),
        "block_n50_variants": n50(lengths_variants),
        "total_bp_phased": sum(lengths_bp),
        "heterozygous_variants": all_row.get("heterozygous_variants", "NA"),
        "heterozygous_snvs": all_row.get("heterozygous_snvs", "NA"),
        "phased_snvs": all_row.get("phased_snvs", "NA"),
        "phased_fraction": all_row.get("phased_fraction", "NA"),
        "phased_snvs_fraction": all_row.get("phased_snvs_fraction", "NA"),
    }

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
