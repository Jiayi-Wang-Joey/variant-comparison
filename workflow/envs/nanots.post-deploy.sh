#!/usr/bin/env bash
set -euo pipefail

# fix: region blocks overlap by one base, duplicating boundary SNVs
d=$(python -c "import nanoTS, os; print(os.path.dirname(nanoTS.__file__))")
sed -i 's/end = min(start + block_size, end_pos)/end = min(start + block_size - 1, end_pos)/' "$d/bam2tpileup.py"
grep -q 'end = min(start + block_size - 1, end_pos)' "$d/bam2tpileup.py"

# fix: count =/X CIGAR ops (pbmm2) as matches
sed -i 's/for op, length in read.cigartuples if op == 0)/for op, length in read.cigartuples if op in (0, 7, 8))/' "$d/snp_feature_extraction_multiple.py"
grep -q 'if op in (0, 7, 8))' "$d/snp_feature_extraction_multiple.py"
