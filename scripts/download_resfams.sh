#!/usr/bin/env bash
# Download and press the Resfams HMM database for IntegronFinder functional annotation.
# Usage: scripts/download_resfams.sh [output_dir]   (requires wget and hmmer)
set -euo pipefail
OUT="${1:-databases/resfams}"
mkdir -p "$OUT"
wget -q -O "$OUT/Resfams.hmm.gz" http://dantaslab.wustl.edu/resfams/Resfams.hmm.gz
gunzip -f "$OUT/Resfams.hmm.gz"
hmmpress -f "$OUT/Resfams.hmm"
echo "Resfams ready in $OUT  ->  use --func_hmm_dir $OUT"
