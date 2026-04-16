#!/usr/bin/env bash
# Run ultra-low VAF detection limit experiments for SV types.
#
# For each config in the ultra-low VAF directory:
#   1. Simulate reads with varforge.
#   2. Run kam on the simulated data.
#   3. Score against truth VCFs.
#
# Outputs sensitivity per VAF per SV type at the end.
#
# Usage:
#   sv_detection_limits.sh --configs-dir DIR --output-dir DIR --reference REF [--force]
#
# Arguments:
#   --configs-dir   Directory containing ultra-low VAF varforge configs.
#   --output-dir    Root output directory for simulation and kam results.
#   --reference     Path to the reference FASTA.
#   --force         Re-run even if output files already exist.
#   -h, --help      Show this help message.
#
# Requires: varforge, kam binary in target/release/kam, Python 3.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
KAM="${REPO}/target/release/kam"

DATA_DIR="${REPO}/docs/benchmarking/03-sv-extended/data"
TARGETS="${DATA_DIR}/sv_suite_targets.fa"

CONFIGS_DIR=""
OUTPUT_DIR=""
REFERENCE=""
FORCE=false

usage() {
    sed -n '2,/^$/s/^# \?//p' "$0"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --configs-dir) CONFIGS_DIR="$2"; shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
        --reference)   REFERENCE="$2"; shift 2 ;;
        --force)       FORCE=true; shift ;;
        -h|--help)     usage 0 ;;
        *)             echo "[ERROR] Unknown argument: $1" >&2; usage 1 ;;
    esac
done

if [[ -z "$CONFIGS_DIR" || -z "$OUTPUT_DIR" || -z "$REFERENCE" ]]; then
    echo "[ERROR] --configs-dir, --output-dir, and --reference are required." >&2
    usage 1
fi

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    exit 1
fi

if ! command -v varforge &>/dev/null; then
    echo "[ERROR] varforge not found in PATH." >&2
    exit 1
fi

if [[ ! -d "$CONFIGS_DIR" ]]; then
    echo "[ERROR] Configs directory not found: $CONFIGS_DIR" >&2
    exit 1
fi

FAILED=()

# Collect all YAML configs.
CONFIGS=()
while IFS= read -r -d '' f; do
    CONFIGS+=("$f")
done < <(find "$CONFIGS_DIR" -name "*.yaml" -print0 | sort -z)

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "[ERROR] No YAML configs found in $CONFIGS_DIR" >&2
    exit 1
fi

echo "Found ${#CONFIGS[@]} configs in ${CONFIGS_DIR}"
echo ""

SIM_DIR="${OUTPUT_DIR}/sim"
KAM_DIR="${OUTPUT_DIR}/kam"
SCORE_DIR="${OUTPUT_DIR}/scores"
mkdir -p "$SIM_DIR" "$KAM_DIR" "$SCORE_DIR"

POSITION_TOLERANCE=10

for cfg in "${CONFIGS[@]}"; do
    basename_cfg=$(basename "$cfg" .yaml)

    echo "=== Processing: ${basename_cfg} ==="

    # Step 1: Simulate with varforge.
    local_sim_dir="${SIM_DIR}/${basename_cfg}"
    if [[ "$FORCE" == "false" ]] && ls "${local_sim_dir}"/*_R1.fastq.gz &>/dev/null 2>&1; then
        echo "  [SKIP] Simulation already exists."
    else
        echo "  [SIM] Running varforge..."
        if ! varforge simulate -c "$cfg" 2>&1 | tail -3; then
            echo "  [ERROR] varforge failed for ${basename_cfg}" >&2
            FAILED+=("sim_${basename_cfg}")
            continue
        fi
    fi

    # Step 2: Run kam.
    local_kam_dir="${KAM_DIR}/${basename_cfg}"
    local_vcf="${local_kam_dir}/variants.vcf"
    if [[ "$FORCE" == "false" && -f "$local_vcf" ]]; then
        echo "  [SKIP] kam output already exists."
    else
        r1=$(ls "${local_sim_dir}"/*_R1.fastq.gz 2>/dev/null | head -1)
        r2=$(ls "${local_sim_dir}"/*_R2.fastq.gz 2>/dev/null | head -1)
        if [[ -z "$r1" || -z "$r2" ]]; then
            echo "  [ERROR] No FASTQs in ${local_sim_dir}" >&2
            FAILED+=("kam_${basename_cfg}")
            continue
        fi

        mkdir -p "$local_kam_dir"
        echo "  [KAM] Running kam..."
        if ! "$KAM" run \
                --r1 "$r1" --r2 "$r2" \
                --targets "$TARGETS" \
                --output-dir "$local_kam_dir" \
                --output-format vcf,tsv \
                2>&1 | grep -E "^\[run\]|ERROR"; then
            echo "  [ERROR] kam failed for ${basename_cfg}" >&2
            FAILED+=("kam_${basename_cfg}")
            continue
        fi
    fi

    # Step 3: Score against truth VCF.
    truth_vcf=$(ls "${local_sim_dir}"/*.truth.vcf 2>/dev/null | head -1)
    if [[ -z "$truth_vcf" ]]; then
        echo "  [WARN] No truth VCF found in ${local_sim_dir}" >&2
        continue
    fi

    # Simple position-based scoring inline.
    tp=0
    fn=0
    fp=0

    # Extract truth positions.
    truth_positions=()
    while IFS=$'\t' read -r _ pos _rest; do
        [[ "$pos" =~ ^[0-9]+$ ]] && truth_positions+=("$pos")
    done < <(grep -v "^#" "$truth_vcf" 2>/dev/null || true)

    # Extract PASS call positions.
    called_positions=()
    if [[ -f "$local_vcf" ]]; then
        while IFS=$'\t' read -r _ pos _ _ _ _ filt _rest; do
            if [[ "$filt" == "PASS" && "$pos" =~ ^[0-9]+$ ]]; then
                called_positions+=("$pos")
            fi
        done < <(grep -v "^#" "$local_vcf" 2>/dev/null || true)
    fi

    # Score: TP if any called position is within tolerance of truth position.
    for t_pos in "${truth_positions[@]}"; do
        matched=false
        for c_pos in "${called_positions[@]}"; do
            diff=$(( c_pos - t_pos ))
            [[ $diff -lt 0 ]] && diff=$(( -diff ))
            if [[ $diff -le $POSITION_TOLERANCE ]]; then
                matched=true
                break
            fi
        done
        if $matched; then
            tp=$((tp + 1))
        else
            fn=$((fn + 1))
        fi
    done

    # FP: called positions not near any truth position.
    for c_pos in "${called_positions[@]}"; do
        near_truth=false
        for t_pos in "${truth_positions[@]}"; do
            diff=$(( c_pos - t_pos ))
            [[ $diff -lt 0 ]] && diff=$(( -diff ))
            if [[ $diff -le $POSITION_TOLERANCE ]]; then
                near_truth=true
                break
            fi
        done
        if ! $near_truth; then
            fp=$((fp + 1))
        fi
    done

    total=$((tp + fn))
    if [[ $total -gt 0 ]]; then
        # Bash integer arithmetic for percentage display.
        sens_pct=$((tp * 100 / total))
    else
        sens_pct=0
    fi

    echo "  [SCORE] TP=${tp} FP=${fp} FN=${fn} Sensitivity=${sens_pct}%"
    echo -e "${basename_cfg}\t${tp}\t${fp}\t${fn}\t${sens_pct}" >> "${SCORE_DIR}/detection_limits.tsv"
done

# Final summary.
echo ""
echo "=== DETECTION LIMITS SUMMARY ==="
echo ""
if [[ -f "${SCORE_DIR}/detection_limits.tsv" ]]; then
    printf "%-40s  %4s  %4s  %4s  %6s\n" "SAMPLE" "TP" "FP" "FN" "SENS%"
    printf "%-40s  %4s  %4s  %4s  %6s\n" "------" "----" "----" "----" "------"
    while IFS=$'\t' read -r name tp fp fn sens; do
        printf "%-40s  %4s  %4s  %4s  %5s%%\n" "$name" "$tp" "$fp" "$fn" "$sens"
    done < "${SCORE_DIR}/detection_limits.tsv"
fi

echo ""
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
echo ""
echo "Detailed scores: ${SCORE_DIR}/detection_limits.tsv"
