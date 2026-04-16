#!/usr/bin/env bash
# Sweep k-mer sizes for SV benchmark datasets and collect results.
#
# Runs kam with each k-mer size on all SV benchmark samples, collecting
# output into per-k subdirectories. Prints a summary table at the end.
#
# Usage:
#   sv_kmer_sweep.sh --data-dir DIR --output-dir DIR --reference REF [--force]
#
# Arguments:
#   --data-dir      Directory containing varforge-simulated SV data.
#                   Expected structure: sim_{type}_vaf{tag}_{rep}/ subdirectories.
#   --output-dir    Root output directory. Results go into kmer_sweep/k{K}/.
#   --reference     Path to the reference FASTA used for targets.
#   --force         Re-run even if output files already exist.
#   -h, --help      Show this help message.
#
# Requires: kam binary in target/release/kam.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
KAM="${REPO}/target/release/kam"

KMER_SIZES=(21 25 27 31 35 41)
SV_TYPES=(sv ins invdel)
VAF_TAGS=(0005 0010 0015 0020 0025 0030 0035 0040 0050 0060 0075 0100 0125 0150 0175 0200 0250 0300 0350 0400 0500 0600 0700 0800 1000)

DATA_DIR=""
OUTPUT_DIR=""
REFERENCE=""
FORCE=false

usage() {
    sed -n '2,/^$/s/^# \?//p' "$0"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)   DATA_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --reference)  REFERENCE="$2"; shift 2 ;;
        --force)      FORCE=true; shift ;;
        -h|--help)    usage 0 ;;
        *)            echo "[ERROR] Unknown argument: $1" >&2; usage 1 ;;
    esac
done

if [[ -z "$DATA_DIR" || -z "$OUTPUT_DIR" || -z "$REFERENCE" ]]; then
    echo "[ERROR] --data-dir, --output-dir, and --reference are required." >&2
    usage 1
fi

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    exit 1
fi

if [[ ! -d "$DATA_DIR" ]]; then
    echo "[ERROR] Data directory not found: $DATA_DIR" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "[ERROR] Reference file not found: $REFERENCE" >&2
    exit 1
fi

# Track failures for final report.
FAILED=()
declare -A RESULTS  # key: "k:type:tag:rep" → "tp:fp:fn"

run_one() {
    local k="$1" sv_type="$2" tag="$3" rep="$4"

    local sim_dir="${DATA_DIR}/sim_${sv_type}_vaf${tag}_${rep}"
    local out_dir="${OUTPUT_DIR}/kmer_sweep/k${k}/${sv_type}_vaf${tag}_${rep}"

    local result_vcf="${out_dir}/variants.vcf"
    if [[ "$FORCE" == "false" && -f "$result_vcf" ]]; then
        return 0
    fi

    local r1 r2
    r1=$(ls "${sim_dir}"/*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sim_dir}"/*_R2.fastq.gz 2>/dev/null | head -1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        return 0  # skip silently if data not yet simulated
    fi

    mkdir -p "$out_dir"

    echo "[RUN] k=${k} ${sv_type} vaf${tag}_${rep}"
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$REFERENCE" \
            -k "$k" \
            --output-dir "$out_dir" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] Failed: k=${k} ${sv_type}_${tag}_${rep}" >&2
        FAILED+=("k${k}_${sv_type}_${tag}_${rep}")
    fi
}

# Main sweep loop.
for k in "${KMER_SIZES[@]}"; do
    echo ""
    echo "========== k-mer size: ${k} =========="
    for sv_type in "${SV_TYPES[@]}"; do
        for tag in "${VAF_TAGS[@]}"; do
            for rep in a b; do
                run_one "$k" "$sv_type" "$tag" "$rep"
            done
        done
    done
done

# Summary table.
echo ""
echo "=== K-MER SWEEP SUMMARY ==="
echo ""
printf "%-6s  %-8s  %-6s  %s\n" "K" "TYPE" "TAG" "STATUS"
printf "%-6s  %-8s  %-6s  %s\n" "---" "--------" "------" "------"

for k in "${KMER_SIZES[@]}"; do
    for sv_type in "${SV_TYPES[@]}"; do
        found=0
        missing=0
        for tag in "${VAF_TAGS[@]}"; do
            for rep in a b; do
                vcf="${OUTPUT_DIR}/kmer_sweep/k${k}/${sv_type}_vaf${tag}_${rep}/variants.vcf"
                if [[ -f "$vcf" ]]; then
                    found=$((found + 1))
                else
                    missing=$((missing + 1))
                fi
            done
        done
        total=$((found + missing))
        printf "%-6s  %-8s  %d/%d complete\n" "$k" "$sv_type" "$found" "$total"
    done
done

echo ""
echo "Failed runs: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
echo ""
echo "Results in: ${OUTPUT_DIR}/kmer_sweep/"
