#!/usr/bin/env bash
# Sweep over read depths from 200K to 10M to measure sensitivity vs depth.
#
# Runs the full titration batch at each depth and writes results to
# benchmarking/results/read_depth/<Nreads>/titration_results_<N>reads.tsv
#
# Usage:
#   bash benchmarking/scripts/sweep_read_depth.sh
#
# Environment:
#   KAM_FASTQ_DIR  — path to TWIST_STDV2_* FASTQ pairs (required if not default)
#
# Runtime estimate (k=31, 100bp targets, single core):
#   200K reads:  ~3 min
#   1M reads:    ~12 min
#   2M reads:    ~25 min  (baseline)
#   5M reads:    ~65 min
#   10M reads:   ~130 min (Phase 2 clustering makes this feasible)

set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SCRIPT="$REPO/benchmarking/scripts/run_titration_batch.py"
RESULTS_BASE="$REPO/benchmarking/results/read_depth"
TARGETS="$REPO/benchmarking/scripts/targets_100bp.fa"
TRUTH_VCF="$REPO/benchmarking/scripts/truth_variants.vcf"
KAM="${KAM_BIN:-$REPO/target/release/kam}"
FASTQ_DIR="${KAM_FASTQ_DIR:-/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs}"

# Read depths to sweep (read pairs per sample).
DEPTHS=(200000 500000 1000000 2000000 5000000 10000000)

echo "=== Read depth sweep ==="
echo "Binary:   $KAM"
echo "Targets:  $TARGETS"
echo "FASTQ:    $FASTQ_DIR"
echo "Results:  $RESULTS_BASE"
echo ""

for depth in "${DEPTHS[@]}"; do
    if (( depth >= 1000000 )); then
        label="$(( depth / 1000000 ))m"
    else
        label="$(( depth / 1000 ))k"
    fi

    echo "--- Depth: $label reads ---"
    results_dir="$RESULTS_BASE/${label}reads"
    mkdir -p "$results_dir"

    # Raise RSS limit to 12GB for deep runs (10M reads ~4GB assembler memory).
    rss_gb=8
    if (( depth >= 5000000 )); then
        rss_gb=12
    fi

    python3 "$SCRIPT" \
        --fastq-dir "$FASTQ_DIR" \
        --targets "$TARGETS" \
        --truth-vcf "$TRUTH_VCF" \
        --kam-binary "$KAM" \
        --results-dir "$results_dir" \
        --reads "$depth" \
        --rss-limit-gb "$rss_gb" \
        --max-vaf 0.35 \
        --min-family-size 2 \
        --target-variants "$TRUTH_VCF" \
        2>&1 | tee "$results_dir/sweep.log"

    echo "Done: $label — results in $results_dir"
    echo ""
done

echo "=== Sweep complete ==="
echo ""
echo "Summary (15ng 2% VAF sensitivity per depth):"
for depth in "${DEPTHS[@]}"; do
    if (( depth >= 1000000 )); then
        label="$(( depth / 1000000 ))m"
    else
        label="$(( depth / 1000 ))k"
    fi
    tsv="$RESULTS_BASE/${label}reads/titration_results_${label}reads.tsv"
    if [[ -f "$tsv" ]]; then
        awk -F'\t' -v d="$label" \
            'NR>1 && $2=="15ng" && $3=="2.0" {printf "  %6s reads: sens=%.3f  SNV=%.3f  indel=%.3f\n", d, $12, $22, $29}' \
            "$tsv"
    fi
done
