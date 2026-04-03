#!/usr/bin/env bash
# Run kam on all 24 titration samples in tumour-informed (monitoring) mode.
#
# Expects the titration FASTQ files in:
#   /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs/
#
# Outputs per-sample results to:
#   docs/benchmarking/04-comparison/kam_results/{sample_id}/
#
# Sample IDs match the format in alignment_baseline.csv:
#   TWIST_STDV2_{ng}ng_VAF_{vaf}pc_DEDUPED_70bp-targets_results
#
# The kam binary is built from the workspace root. Build it first:
#   cargo build --release
#
# Usage:
#   bash run_kam_titration.sh [--dry-run]
#
# Options:
#   --dry-run   Print commands without executing them.
#
# Environment overrides:
#   KAM_FASTQ_DIR   Path to FASTQ files (default: /mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs)
#   KAM_BINARY      Path to kam binary (default: <repo>/target/release/kam)

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        *) echo "[ERROR] Unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

KAM="${KAM_BINARY:-${REPO}/target/release/kam}"
FASTQ_DIR="${KAM_FASTQ_DIR:-/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs}"
TARGETS="${SCRIPT_DIR}/../../../benchmarking/snvindel/scripts/targets_100bp.fa"
TRUTH_VCF="${SCRIPT_DIR}/../../../benchmarking/snvindel/scripts/truth_variants.vcf"
KAM_RESULTS_DIR="${SCRIPT_DIR}/../kam_results"
PERF_TSV="${KAM_RESULTS_DIR}/perf.tsv"

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if [[ $DRY_RUN -eq 0 ]]; then
    if [[ ! -x "$KAM" ]]; then
        echo "[ERROR] kam binary not found or not executable: $KAM" >&2
        echo "        Build with: cargo build --release" >&2
        exit 1
    fi

    if [[ ! -f "$TARGETS" ]]; then
        echo "[ERROR] Targets FASTA not found: $TARGETS" >&2
        exit 1
    fi

    if [[ ! -f "$TRUTH_VCF" ]]; then
        echo "[ERROR] Truth VCF not found: $TRUTH_VCF" >&2
        exit 1
    fi

    if [[ ! -d "$FASTQ_DIR" ]]; then
        echo "[ERROR] FASTQ directory not found: $FASTQ_DIR" >&2
        exit 1
    fi
fi

mkdir -p "$KAM_RESULTS_DIR"

# Write performance log header if it does not exist yet.
if [[ ! -f "$PERF_TSV" ]]; then
    printf "sample_id\twall_seconds\tmax_rss_kb\n" > "$PERF_TSV"
fi

# ---------------------------------------------------------------------------
# Helper: derive sample_id from FASTQ filename
#
# Input:  TWIST_STDV2_15ng_VAF_0p1pc_22KVL2LT3_CCAACTCCGA-GGAGTAACGC_L008_R1.fastq.gz
# Output: TWIST_STDV2_15ng_VAF_0p1pc_DEDUPED_70bp-targets_results
#
# The sample_id mirrors the format used in alignment_baseline.csv so that the
# concordance script can join the two tables without a separate mapping file.
# ---------------------------------------------------------------------------
derive_sample_id() {
    local filename="$1"
    local base
    base="$(basename "$filename")"

    # Extract fields 3 (ng) and 5 (vaf+pc) from underscore-split filename.
    # TWIST_STDV2_{ng}ng_VAF_{vaf}pc_{instrument}_{barcodes}_{lane}_R1.fastq.gz
    #   idx: 1     2  3    4  5    6         7        8       9 10
    local ng_field vaf_field
    ng_field="$(echo "$base" | cut -d_ -f3)"    # e.g. 15ng
    vaf_field="$(echo "$base" | cut -d_ -f5)"   # e.g. 0p1pc

    echo "TWIST_STDV2_${ng_field}_VAF_${vaf_field}_DEDUPED_70bp-targets_results"
}

# ---------------------------------------------------------------------------
# Process each R1 file
# ---------------------------------------------------------------------------
shopt -s nullglob
r1_files=("${FASTQ_DIR}"/TWIST_STDV2_*_R1.fastq.gz)

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "[ERROR] No R1 FASTQ files found in ${FASTQ_DIR}" >&2
    exit 1
fi

echo "[INFO] Found ${#r1_files[@]} R1 files." >&2

for r1 in "${r1_files[@]}"; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"

    if [[ ! -f "$r2" ]]; then
        echo "[WARN] R2 not found for $(basename "$r1"); skipping." >&2
        continue
    fi

    sample_id="$(derive_sample_id "$r1")"
    out_dir="${KAM_RESULTS_DIR}/${sample_id}"

    echo "[INFO] Processing ${sample_id}..." >&2

    cmd=(
        "$KAM" run
            --r1 "$r1"
            --r2 "$r2"
            --targets "$TARGETS"
            --output-dir "$out_dir"
            --output-format vcf,tsv
            --target-variants "$TRUTH_VCF"
    )

    if [[ $DRY_RUN -eq 1 ]]; then
        echo "  DRY-RUN: ${cmd[*]}" >&2
        continue
    fi

    mkdir -p "$out_dir"
    time_log="$(mktemp)"

    /usr/bin/time -v "${cmd[@]}" 2>"$time_log" || {
        echo "[ERROR] kam failed for ${sample_id}" >&2
        cat "$time_log" >&2
        rm -f "$time_log"
        continue
    }

    # Parse /usr/bin/time -v output for wall time and peak RSS.
    wall_sec="$(grep -oP '(?<=Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): )\S+' "$time_log" \
        | awk -F: '{
            if (NF==3) { print $1*3600 + $2*60 + $3 }
            else       { print $1*60   + $2 }
          }')"
    max_rss="$(grep -oP '(?<=Maximum resident set size \(kbytes\): )\d+' "$time_log")"
    rm -f "$time_log"

    printf "%s\t%s\t%s\n" \
        "$sample_id" "${wall_sec:-NA}" "${max_rss:-NA}" \
        >> "$PERF_TSV"

    echo "[INFO] Done: ${sample_id} (${wall_sec:-?}s, ${max_rss:-?} kB RSS)" >&2
done

echo "[INFO] All samples processed. Results in: ${KAM_RESULTS_DIR}" >&2
