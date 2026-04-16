#!/usr/bin/env bash
# Sweep confidence thresholds and min-alt-molecule counts for SV calling.
#
# Indexes once per sample, then calls with each parameter combination.
# This avoids redundant indexing and speeds up the sweep significantly.
#
# Usage:
#   sv_threshold_sweep.sh --data-dir DIR --output-dir DIR [--force]
#
# Arguments:
#   --data-dir      Directory containing varforge-simulated SV data.
#                   Expected structure: sim_{type}_vaf{tag}_{rep}/ subdirectories.
#   --output-dir    Root output directory. Results go into
#                   threshold_sweep/conf{CONF}_alt{ALT}/.
#   --force         Re-run even if output files already exist.
#   -h, --help      Show this help message.
#
# Requires: kam binary in target/release/kam.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
KAM="${REPO}/target/release/kam"

CONFIDENCE_LEVELS=(0.80 0.90 0.95 0.99)
ALT_MOLECULE_COUNTS=(1 2 3)
SV_TYPES=(sv ins invdel)
VAF_TAGS=(0005 0010 0015 0020 0025 0030 0035 0040 0050 0060 0075 0100 0125 0150 0175 0200 0250 0300 0350 0400 0500 0600 0700 0800 1000)

DATA_DIR=""
OUTPUT_DIR=""
FORCE=false

usage() {
    sed -n '2,/^$/s/^# \?//p' "$0"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)   DATA_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --force)      FORCE=true; shift ;;
        -h|--help)    usage 0 ;;
        *)            echo "[ERROR] Unknown argument: $1" >&2; usage 1 ;;
    esac
done

if [[ -z "$DATA_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "[ERROR] --data-dir and --output-dir are required." >&2
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

FAILED=()

# Format confidence for directory name: 0.95 → 095.
fmt_conf() {
    local c="$1"
    printf "%03d" "$(echo "$c * 100" | bc | cut -d. -f1)"
}

call_one() {
    local sv_type="$1" tag="$2" rep="$3" conf="$4" alt="$5"

    local conf_tag
    conf_tag=$(fmt_conf "$conf")
    local sweep_dir="${OUTPUT_DIR}/threshold_sweep/conf${conf_tag}_alt${alt}"
    local out_dir="${sweep_dir}/${sv_type}_vaf${tag}_${rep}"
    local result_vcf="${out_dir}/variants.vcf"

    if [[ "$FORCE" == "false" && -f "$result_vcf" ]]; then
        return 0
    fi

    local sim_dir="${DATA_DIR}/sim_${sv_type}_vaf${tag}_${rep}"
    local r1 r2
    r1=$(ls "${sim_dir}"/*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sim_dir}"/*_R2.fastq.gz 2>/dev/null | head -1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        return 0
    fi

    mkdir -p "$out_dir"

    echo "[CALL] conf=${conf} alt=${alt} ${sv_type} vaf${tag}_${rep}"
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "${DATA_DIR}/../data/sv_suite_targets.fa" \
            --sv-min-confidence "$conf" \
            --sv-min-alt-molecules "$alt" \
            --output-dir "$out_dir" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] Failed: conf=${conf} alt=${alt} ${sv_type}_${tag}_${rep}" >&2
        FAILED+=("conf${conf_tag}_alt${alt}_${sv_type}_${tag}_${rep}")
    fi
}

# Main sweep loop.
for conf in "${CONFIDENCE_LEVELS[@]}"; do
    for alt in "${ALT_MOLECULE_COUNTS[@]}"; do
        conf_tag=$(fmt_conf "$conf")
        echo ""
        echo "========== conf=${conf} alt=${alt} =========="
        for sv_type in "${SV_TYPES[@]}"; do
            for tag in "${VAF_TAGS[@]}"; do
                for rep in a b; do
                    call_one "$sv_type" "$tag" "$rep" "$conf" "$alt"
                done
            done
        done
    done
done

# Summary table.
echo ""
echo "=== THRESHOLD SWEEP SUMMARY ==="
echo ""
printf "%-6s  %-4s  %-8s  %s\n" "CONF" "ALT" "TYPE" "STATUS"
printf "%-6s  %-4s  %-8s  %s\n" "------" "----" "--------" "------"

for conf in "${CONFIDENCE_LEVELS[@]}"; do
    conf_tag=$(fmt_conf "$conf")
    for alt in "${ALT_MOLECULE_COUNTS[@]}"; do
        for sv_type in "${SV_TYPES[@]}"; do
            found=0
            missing=0
            for tag in "${VAF_TAGS[@]}"; do
                for rep in a b; do
                    vcf="${OUTPUT_DIR}/threshold_sweep/conf${conf_tag}_alt${alt}/${sv_type}_vaf${tag}_${rep}/variants.vcf"
                    if [[ -f "$vcf" ]]; then
                        found=$((found + 1))
                    else
                        missing=$((missing + 1))
                    fi
                done
            done
            total=$((found + missing))
            printf "%-6s  %-4s  %-8s  %d/%d complete\n" "$conf" "$alt" "$sv_type" "$found" "$total"
        done
    done
done

echo ""
echo "Failed runs: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
echo ""
echo "Results in: ${OUTPUT_DIR}/threshold_sweep/"
