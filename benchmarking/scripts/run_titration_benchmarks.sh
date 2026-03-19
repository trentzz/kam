#!/usr/bin/env bash
# Run kam on all titration FASTQ pairs and score against truth variants.
#
# Filename pattern:
#   TWIST_STDV2_{ng}ng_VAF_{vaf}pc_{instrument}_{barcodes}_{lane}_R{1,2}.fastq.gz
#
# Sample names in all output are anonymized to:
#   Sample_{ng}ng_VAF_{vaf}pc
# No instrument ID, barcode, or lane information appears in output files.
#
# Results are appended to $RESULTS_DIR/titration_results.tsv.
# Timing and memory go to $RESULTS_DIR/titration_perf.tsv.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
KAM=/home/trent/tfiles/code/kam/target/release/kam
TARGETS=/home/trent/tfiles/code/kam/benchmarking/scripts/targets_100bp.fa
FASTQ_DIR=/mnt/tzeng-local/tzeng-thesis/titration-nondedup/fastqs
RESULTS_DIR=/home/trent/tfiles/code/kam/benchmarking/results/tables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCORE_PY="${SCRIPT_DIR}/score_variants.py"
TRUTH_VCF="${SCRIPT_DIR}/truth_variants.vcf"

# Output subdirectory for kam runs (under RESULTS_DIR, never contains raw IDs).
KAM_RUNS_DIR="${RESULTS_DIR}/titration_runs"

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found or not executable: $KAM" >&2
    exit 1
fi

if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Targets FASTA not found: $TARGETS" >&2
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "[ERROR] FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

if [[ ! -f "$TRUTH_VCF" ]]; then
    echo "[INFO] truth_variants.vcf not found; generating from TSV..." >&2
    python3 "${SCRIPT_DIR}/tsv_to_vcf.py" \
        --input /mnt/tzeng-local/tzeng-thesis/titration.probes.QC.pass.tsv \
        --output "$TRUTH_VCF"
fi

mkdir -p "$RESULTS_DIR" "$KAM_RUNS_DIR"

RESULTS_TSV="${RESULTS_DIR}/titration_results.tsv"
PERF_TSV="${RESULTS_DIR}/titration_perf.tsv"

if [[ ! -f "$PERF_TSV" ]]; then
    printf "sample_name\twall_seconds\tmax_rss_kb\n" > "$PERF_TSV"
fi

# ---------------------------------------------------------------------------
# Process each R1 file
# ---------------------------------------------------------------------------
# Pattern: TWIST_STDV2_{ng}ng_VAF_{vaf}pc_{rest}_R1.fastq.gz
#   vaf may use 'p' as a decimal separator: 0p1 = 0.1, 0p001 = 0.001
# ---------------------------------------------------------------------------

shopt -s nullglob
r1_files=("${FASTQ_DIR}"/TWIST_STDV2_*_R1.fastq.gz)

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "[ERROR] No R1 FASTQ files found in ${FASTQ_DIR}" >&2
    exit 1
fi

echo "[INFO] Found ${#r1_files[@]} R1 files in ${FASTQ_DIR}" >&2

for r1 in "${r1_files[@]}"; do
    r1_base="$(basename "$r1")"

    # Derive R2 path.
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    if [[ ! -f "$r2" ]]; then
        echo "[WARN] R2 not found for ${r1_base}; skipping." >&2
        continue
    fi

    # Extract ng amount: second field after splitting on '_'.
    # Filename: TWIST_STDV2_{ng}ng_VAF_{vaf}pc_{rest}
    #           field:  0     1   2    3  4    5...
    ng_field="$(echo "$r1_base" | cut -d_ -f3)"          # e.g. 15ng
    ng="${ng_field%ng}"                                    # e.g. 15

    # VAF field: TWIST_STDV2_{ng}ng_VAF_{vaf}pc_{rest}
    #            field 5 is the vaf+pc string.
    vaf_raw="$(echo "$r1_base" | cut -d_ -f5)"            # e.g. 0p1pc
    vaf_stripped="${vaf_raw%pc}"                           # e.g. 0p1
    vaf_numeric="${vaf_stripped/p/.}"                      # e.g. 0.1

    # Anonymized sample name: no instrument, barcodes, or lane.
    sample_name="Sample_${ng}ng_VAF_${vaf_stripped}pc"

    echo "[INFO] Processing ${sample_name}..." >&2

    kam_out="${KAM_RUNS_DIR}/${sample_name}"
    mkdir -p "$kam_out"

    time_log="$(mktemp)"

    /usr/bin/time -v \
        "$KAM" run \
            --r1 "$r1" \
            --r2 "$r2" \
            --targets "$TARGETS" \
            --output-dir "$kam_out" \
            --output-format vcf,tsv \
        2>"$time_log" || {
            echo "[ERROR] kam failed for ${sample_name}" >&2
            cat "$time_log" >&2
            rm -f "$time_log"
            continue
        }

    # Parse /usr/bin/time -v output.
    wall_sec="$(grep -oP '(?<=Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): )\S+' "$time_log" \
        | awk -F: '{
            if (NF==3) { print $1*3600 + $2*60 + $3 }
            else       { print $1*60   + $2 }
          }')"
    max_rss="$(grep -oP '(?<=Maximum resident set size \(kbytes\): )\d+' "$time_log")"
    rm -f "$time_log"

    printf "%s\t%s\t%s\n" \
        "$sample_name" "${wall_sec:-NA}" "${max_rss:-NA}" \
        >> "$PERF_TSV"

    # Locate kam VCF output.
    called_vcf="$(find "$kam_out" -maxdepth 2 -name "*.vcf" 2>/dev/null | sort | head -1)"

    if [[ -z "$called_vcf" ]]; then
        echo "[WARN] No VCF output found for ${sample_name}; skipping scoring." >&2
        continue
    fi

    python3 "$SCORE_PY" \
        --truth "$TRUTH_VCF" \
        --called "$called_vcf" \
        --output "$RESULTS_TSV" \
        --sample-name "$sample_name" \
        --vaf "$vaf_numeric" \
        --dna-input "$ng"
done

echo "[INFO] Titration benchmarks complete. Results: ${RESULTS_TSV}" >&2
