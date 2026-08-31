#!/usr/bin/env python3
# convert_isolaser_vcf.py
# Usage: python3 convert_isolaser_vcf.py input.vcf.gz sample_name output.vcf.gz

import pysam
import sys

input_vcf  = sys.argv[1]
sample_name = sys.argv[2]
output_vcf = sys.argv[3]

vcf_in = pysam.VariantFile(input_vcf)

# Build new header with PS tag and fixed sample name
new_header = pysam.VariantHeader()
for rec in vcf_in.header.records:
    new_header.add_record(rec)
if "PS" not in vcf_in.header.formats:
    new_header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set Identifier">')
new_header.add_sample(sample_name)

# First pass: deduplicate by position
# Priority: phased record > unphased; then highest QUAL within same phasing status
def is_phased(rec):
    old_sample = list(vcf_in.header.samples)[0]
    return rec.samples[old_sample].phased

records_by_pos = {}
for rec in vcf_in:
    key = (rec.chrom, rec.pos)
    if key not in records_by_pos:
        records_by_pos[key] = rec
    else:
        existing = records_by_pos[key]
        rec_phased = is_phased(rec)
        existing_phased = is_phased(existing)
        # prefer phased over unphased
        if rec_phased and not existing_phased:
            records_by_pos[key] = rec
        # if same phasing status, keep higher QUAL
        elif rec_phased == existing_phased:
            if rec.qual is not None and (existing.qual is None or rec.qual > existing.qual):
                records_by_pos[key] = rec

# Second pass: assign PS from PG, write output
phase_group_to_ps = {}  # PG group name -> position of first variant

vcf_out = pysam.VariantFile(output_vcf, "w", header=new_header)

for key in sorted(records_by_pos):
    rec = records_by_pos[key]
    old_sample = list(vcf_in.header.samples)[0]

    new_rec = vcf_out.new_record()
    new_rec.chrom  = rec.chrom
    new_rec.pos    = rec.pos
    new_rec.id     = rec.id
    new_rec.alleles = rec.alleles
    new_rec.qual   = rec.qual
    new_rec.filter.add(list(rec.filter)[0])

    # Copy INFO
    for key_info, val in rec.info.items():
        try:
            new_rec.info[key_info] = val
        except Exception:
            pass

    # Copy FORMAT fields + add PS
    s = rec.samples[old_sample]
    ns = new_rec.samples[sample_name]

    for fmt in ["GT", "GQ", "PGT", "DP", "MIN_DP", "AC", "AB", "PL", "PG"]:
        try:
            ns[fmt] = s[fmt]
        except Exception:
            pass

    # Convert PG → PS
    pg = s.get("PG", None)
    if pg is not None:
        parts = pg.split("_")
        # strip trailing haplotype digit (_0 or _1)
        group = "_".join(parts[:-1]) if parts[-1] in ("0", "1") else pg
        if group not in phase_group_to_ps:
            phase_group_to_ps[group] = rec.pos
        ns["PS"] = phase_group_to_ps[group]

    # Preserve phased GT
    ns.phased = s.phased

    vcf_out.write(new_rec)

vcf_in.close()
vcf_out.close()