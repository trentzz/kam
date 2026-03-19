#!/usr/bin/env bash
# Run synthetic benchmarks for kam using varforge-generated data.
#
# For each varforge config in benchmarking/scripts/configs/:
#   - Twist duplex UMI configs  (twist_*.yaml)  are run first.
#   - Other-technology configs  (simplex_umi_panel.yaml, wgs_no_umi.yaml, etc.)
#     are run in a separate section afterwards.
#
# Results are appended to $RESULTS_DIR/synthetic_results.tsv.
# Memory and timing data go to $RESULTS_DIR/synthetic_perf.tsv.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
VARFORGE=/tmp/varforge/target/release/varforge
KAM=/home/trent/tfiles/code/kam/target/release/kam
TARGETS=/home/trent/tfiles/code/kam/benchmarking/scripts/targets_100bp.fa
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/configs"
RESULTS_DIR=/home/trent/tfiles/code/kam/benchmarking/results/tables
SYNTHETIC_DIR=/home/trent/tfiles/code/kam/benchmarking/synthetic_data
SCORE_PY="${SCRIPT_DIR}/score_variants.py"
TRUTH_VCF="${SCRIPT_DIR}/truth_variants.vcf"

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
for bin in "$VARFORGE" "$KAM"; do
    if [[ ! -x "$bin" ]]; then
        echo "[ERROR] Binary not found or not executable: $bin" >&2
        exit 1
    fi
done

if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Targets FASTA not found: $TARGETS" >&2
    exit 1
fi

if [[ ! -f "$TRUTH_VCF" ]]; then
    echo "[INFO] truth_variants.vcf not found; generating from TSV..." >&2
    python3 "${SCRIPT_DIR}/tsv_to_vcf.py" \
        --input /mnt/tzeng-local/tzeng-thesis/titration.probes.QC.pass.tsv \
        --output "$TRUTH_VCF"
fi

mkdir -p "$RESULTS_DIR" "$SYNTHETIC_DIR"

RESULTS_TSV="${RESULTS_DIR}/synthetic_results.tsv"
PERF_TSV="${RESULTS_DIR}/synthetic_perf.tsv"

# Write performance header if file does not exist.
if [[ ! -f "$PERF_TSV" ]]; then
    printf "config\tsample\twall_seconds\tmax_rss_kb\n" > "$PERF_TSV"
fi

# ---------------------------------------------------------------------------
# Helper: run one sample directory
# ---------------------------------------------------------------------------
run_sample() {
    local config_name="$1"   # e.g. twist_sensitivity_sweep
    local sample_dir="$2"    # path to a sample sub-directory produced by varforge
    local section="$3"       # "twist" or "other"

    local sample_name
    sample_name="$(basename "$sample_dir")"

    # Locate R1 and R2 inside the sample directory.
    local r1 r2
    r1="$(find "$sample_dir" -maxdepth 2 -name "*_R1*.fastq.gz" -o -name "*_R1*.fastq" 2>/dev/null | sort | head -1)"
    r2="$(find "$sample_dir" -maxdepth 2 -name "*_R2*.fastq.gz" -o -name "*_R2*.fastq" 2>/dev/null | sort | head -1)"

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[WARN] Could not find R1/R2 in ${sample_dir}; skipping." >&2
        return 0
    fi

    local kam_out="${sample_dir}/kam_out"
    mkdir -p "$kam_out"

    local time_log
    time_log="$(mktemp)"

    echo "[INFO] Running kam on ${config_name}/${sample_name}..." >&2

    /usr/bin/time -v \
        "$KAM" run \
            --r1 "$r1" \
            --r2 "$r2" \
            --targets "$TARGETS" \
            --output-dir "$kam_out" \
            --output-format vcf,tsv \
        2>"$time_log" || {
            echo "[ERROR] kam failed for ${config_name}/${sample_name}" >&2
            cat "$time_log" >&2
            rm -f "$time_log"
            return 1
        }

    # Parse /usr/bin/time -v output.
    local wall_sec max_rss
    wall_sec="$(grep -oP '(?<=Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): )\S+' "$time_log" \
        | awk -F: '{
            if (NF==3) { print $1*3600 + $2*60 + $3 }
            else       { print $1*60   + $2 }
          }')"
    max_rss="$(grep -oP '(?<=Maximum resident set size \(kbytes\): )\d+' "$time_log")"
    rm -f "$time_log"

    printf "%s\t%s\t%s\t%s\n" \
        "$config_name" "$sample_name" "${wall_sec:-NA}" "${max_rss:-NA}" \
        >> "$PERF_TSV"

    # Locate kam VCF output.
    local called_vcf
    called_vcf="$(find "$kam_out" -maxdepth 2 -name "*.vcf" 2>/dev/null | sort | head -1)"

    if [[ -z "$called_vcf" ]]; then
        echo "[WARN] No VCF output found in ${kam_out}; skipping scoring." >&2
        return 0
    fi

    # Extract VAF and DNA input from sample directory name if possible.
    # varforge typically names samples like vaf_0.001 or vaf1pc etc.
    local vaf dna_input
    vaf="$(echo "$sample_name" | grep -oP '(?i)(?:vaf[_-]?)[\d.]+' | grep -oP '[\d.]+' || echo "NA")"
    dna_input="NA"

    python3 "$SCORE_PY" \
        --truth "$TRUTH_VCF" \
        --called "$called_vcf" \
        --output "$RESULTS_TSV" \
        --sample-name "${config_name}_${sample_name}" \
        --vaf "${vaf}" \
        --dna-input "${dna_input}"
}

# ---------------------------------------------------------------------------
# Section 1: Twist duplex UMI configs
# ---------------------------------------------------------------------------
echo "[INFO] === Section 1: Twist duplex UMI configs ===" >&2

for config in "${CONFIG_DIR}"/twist_*.yaml; do
    [[ -f "$config" ]] || continue
    config_name="$(basename "$config" .yaml)"
    output_dir="${SYNTHETIC_DIR}/${config_name}"
    mkdir -p "$output_dir"

    echo "[INFO] Generating synthetic data: ${config_name}..." >&2
    "$VARFORGE" simulate --output-dir "$output_dir" "$config"

    # Each sub-directory is one sample (one VAF level).
    for sample_dir in "${output_dir}"/*/; do
        [[ -d "$sample_dir" ]] || continue
        run_sample "$config_name" "$sample_dir" "twist"
    done
done

# ---------------------------------------------------------------------------
# Section 2: Other-technology configs
# ---------------------------------------------------------------------------
echo "[INFO] === Section 2: Other-technology configs ===" >&2

OTHER_CONFIGS=(
    "${CONFIG_DIR}/simplex_umi_panel.yaml"
    "${CONFIG_DIR}/wgs_no_umi.yaml"
    "${CONFIG_DIR}/cfdna_simplex_panel.yaml"
    "${CONFIG_DIR}/high_depth_amplicon.yaml"
)

for config in "${OTHER_CONFIGS[@]}"; do
    [[ -f "$config" ]] || { echo "[WARN] Config not found: ${config}; skipping." >&2; continue; }
    config_name="$(basename "$config" .yaml)"
    output_dir="${SYNTHETIC_DIR}/${config_name}"
    mkdir -p "$output_dir"

    echo "[INFO] Generating synthetic data: ${config_name}..." >&2
    "$VARFORGE" simulate --output-dir "$output_dir" "$config"

    for sample_dir in "${output_dir}"/*/; do
        [[ -d "$sample_dir" ]] || continue
        run_sample "$config_name" "$sample_dir" "other"
    done
done

echo "[INFO] Synthetic benchmarks complete. Results: ${RESULTS_TSV}" >&2
