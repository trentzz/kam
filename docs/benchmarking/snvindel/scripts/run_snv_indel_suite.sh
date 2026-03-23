#!/usr/bin/env bash
# Run kam on all SNV and indel benchmark datasets.
#
# For each varforge-simulated dataset:
#   1. Discovery mode  — no --target-variants
#   2. Tumour-informed mode — --target-variants using the varforge truth VCF
#
# Outputs per dataset:
#   results/kam_snv_vaf{tag}_{rep}/calls_discovery.vcf
#   results/kam_snv_vaf{tag}_{rep}/calls_tumour_informed.vcf
#   results/kam_indel_vaf{tag}_{rep}/calls_discovery.vcf
#   results/kam_indel_vaf{tag}_{rep}/calls_tumour_informed.vcf
#
# Usage: bash run_snv_indel_suite.sh [--force]
#
# --force: re-run even if output files already exist.
#
# Requires: target/release/kam, Python 3.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
KAM="${REPO}/target/release/kam"
EXPDIR="${REPO}/docs/benchmarking/snvindel"
DATADIR="${EXPDIR}/data"
RESDIR="${EXPDIR}/results"
TARGETS="${DATADIR}/snvindel_targets.fa"
FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    exit 1
fi
if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Targets FASTA not found: $TARGETS" >&2
    exit 1
fi

FAILED=()

run_dataset() {
    local type="$1"   # snv or indel
    local tag="$2"    # e.g. 0100
    local rep="$3"    # a or b

    local sim_dir="${RESDIR}/sim_${type}_vaf${tag}_${rep}"
    local out_dir="${RESDIR}/kam_${type}_vaf${tag}_${rep}"

    local disc_vcf="${out_dir}/calls_discovery.vcf"
    local ti_vcf="${out_dir}/calls_tumour_informed.vcf"

    # Skip if already complete.
    if [[ "$FORCE" == "false" && -f "$disc_vcf" && -f "$ti_vcf" ]]; then
        return 0
    fi

    # Locate R1/R2 FASTQs.
    local r1 r2
    r1=$(ls "${sim_dir}"/*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sim_dir}"/*_R2.fastq.gz 2>/dev/null | head -1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[WARN] No FASTQs in ${sim_dir}" >&2
        FAILED+=("${type}_${tag}_${rep}")
        return 0
    fi

    # Locate varforge truth VCF.
    local truth_vcf
    truth_vcf=$(ls "${sim_dir}"/*.truth.vcf 2>/dev/null | head -1)
    if [[ -z "$truth_vcf" ]]; then
        echo "[WARN] No truth VCF in ${sim_dir}" >&2
        FAILED+=("${type}_${tag}_${rep}")
        return 0
    fi

    mkdir -p "$out_dir"

    echo "[RUN] ${type} vaf${tag}_${rep} — discovery"
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$TARGETS" \
            --output-dir "${out_dir}/tmp_disc" \
            --output-format vcf \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] discovery failed for ${type}_${tag}_${rep}" >&2
        FAILED+=("disc_${type}_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_disc/variants.vcf" "$disc_vcf"
    rm -rf "${out_dir}/tmp_disc"

    echo "[RUN] ${type} vaf${tag}_${rep} — tumour-informed"
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$TARGETS" \
            --target-variants "$truth_vcf" \
            --output-dir "${out_dir}/tmp_ti" \
            --output-format vcf \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] tumour-informed failed for ${type}_${tag}_${rep}" >&2
        FAILED+=("ti_${type}_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_ti/variants.vcf" "$ti_vcf"
    rm -rf "${out_dir}/tmp_ti"
}

VAF_TAGS=(0005 0010 0015 0020 0025 0030 0035 0040 0050 0060 0075 0100 0125 0150 0175 0200 0250 0300 0350 0400 0500 0600 0700 0800 1000)

for tag in "${VAF_TAGS[@]}"; do
    for rep in a b; do
        run_dataset snv   "$tag" "$rep"
        run_dataset indel "$tag" "$rep"
    done
done

echo ""
echo "=== DONE ==="
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
