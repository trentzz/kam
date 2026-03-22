#!/usr/bin/env bash
# Sweep over k-mer sizes from 21 to 41 to measure sensitivity vs k.
#
# Uses the fixed 100bp target panel (101bp actual sequences). All k values
# in this range satisfy the minimum detectable length requirement (2k ≤ 100bp).
# Reads capped at 2M per sample (same as the baseline benchmark).
#
# Results written to benchmarking/results/kmer_size/k<N>/titration_results_2mreads.tsv
#
# Usage:
#   bash benchmarking/scripts/sweep_kmer_size.sh
#
# Environment:
#   KAM_FASTQ_DIR  — path to TWIST_STDV2_* FASTQ pairs (required if not default)
#
# Note: k=31 is the baseline. Smaller k increases spurious branching in the
# de Bruijn graph; larger k reduces it but requires longer reads to cover anchors.
# For 150bp reads, k up to ~60 is geometrically safe, but the maximum useful k
# for 100bp targets is ~40 (2k=80 < 101bp actual sequence length).

set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SCRIPT="$REPO/benchmarking/scripts/run_titration_batch.py"
RESULTS_BASE="$REPO/benchmarking/results/kmer_size"
TARGETS="$REPO/benchmarking/scripts/targets_100bp.fa"
TRUTH_VCF="$REPO/benchmarking/scripts/truth_variants.vcf"
KAM="${KAM_BIN:-$REPO/target/release/kam}"
FASTQ_DIR="${KAM_FASTQ_DIR:-/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs}"
READS=2000000

# K values to sweep. Odd numbers preferred (avoids palindrome edge cases),
# but the canonical function handles even k correctly.
K_VALUES=(21 25 27 29 31 33 35 37 41)

echo "=== K-mer size sweep ==="
echo "Binary:   $KAM"
echo "Targets:  $TARGETS (100bp, 101bp actual sequences)"
echo "FASTQ:    $FASTQ_DIR"
echo "Reads:    $READS per sample"
echo "Results:  $RESULTS_BASE"
echo ""

for k in "${K_VALUES[@]}"; do
    echo "--- k=$k ---"
    results_dir="$RESULTS_BASE/k${k}"
    mkdir -p "$results_dir"

    python3 "$SCRIPT" \
        --fastq-dir "$FASTQ_DIR" \
        --targets "$TARGETS" \
        --truth-vcf "$TRUTH_VCF" \
        --kam-binary "$KAM" \
        --results-dir "$results_dir" \
        --reads "$READS" \
        --max-vaf 0.35 \
        --min-family-size 2 \
        --target-variants "$TRUTH_VCF" \
        --kmer-size "$k" \
        --output "titration_results_k${k}.tsv" \
        2>&1 | tee "$results_dir/sweep.log"

    echo "Done: k=$k — results in $results_dir"
    echo ""
done

echo "=== Sweep complete ==="
echo ""
echo "Summary (15ng 2% VAF sensitivity per k):"
for k in "${K_VALUES[@]}"; do
    tsv="$RESULTS_BASE/k${k}/titration_results_k${k}.tsv"
    if [[ -f "$tsv" ]]; then
        awk -F'\t' -v k="$k" \
            'NR>1 && $2=="15ng" && $3=="2.0" {printf "  k=%-2s: sens=%.3f  SNV=%.3f  indel=%.3f  pathfind_ms=%s\n", k, $12, $22, $29, $35}' \
            "$tsv"
    fi
done
