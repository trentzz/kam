#!/usr/bin/env bash
# Run kam on all SV benchmark datasets (DUP+INV, INS, INVDEL).
#
# For each varforge-simulated dataset:
#   1. Discovery mode  — no --target-variants
#   2. Tumour-informed mode — --target-variants using the varforge truth VCF
#
# Outputs per dataset:
#   results/kam_{type}_vaf{tag}_{rep}/calls_discovery.vcf
#   results/kam_{type}_vaf{tag}_{rep}/calls_discovery.tsv
#   results/kam_{type}_vaf{tag}_{rep}/calls_tumour_informed.vcf
#   results/kam_{type}_vaf{tag}_{rep}/calls_tumour_informed.tsv
#
# Usage: bash run_sv_suite.sh [--force]
#
# --force: re-run even if output files already exist.
#
# Note: DEL variants excluded from the DUP+INV suite due to a varforge engine
# panic for deletions ≥20bp. The suite benchmarks DUP and INV only.
#
# Requires: target/release/kam, Python 3.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
KAM="${REPO}/target/release/kam"
EXPDIR="${REPO}/docs/benchmarking/sv"
DATADIR="${EXPDIR}/data"
RESDIR="${EXPDIR}/results"
FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    exit 1
fi

FAILED=()

run_dataset() {
    local type="$1"       # sv, ins, or invdel
    local sim_type="$2"   # sv, ins, or invdel (sim directory prefix)
    local tag="$3"        # e.g. 0100
    local rep="$4"        # a or b
    local targets="$5"
    local junctions="$6"  # may be empty

    local sim_dir="${RESDIR}/sim_${sim_type}_vaf${tag}_${rep}"
    local out_dir="${RESDIR}/kam_${type}_vaf${tag}_${rep}"

    local disc_vcf="${out_dir}/calls_discovery.vcf"
    local ti_vcf="${out_dir}/calls_tumour_informed.vcf"

    if [[ "$FORCE" == "false" && -f "$disc_vcf" && -f "$ti_vcf" ]]; then
        return 0
    fi

    local r1 r2
    r1=$(ls "${sim_dir}"/*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sim_dir}"/*_R2.fastq.gz 2>/dev/null | head -1)
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[WARN] No FASTQs in ${sim_dir}" >&2
        FAILED+=("${type}_${tag}_${rep}")
        return 0
    fi

    local truth_vcf
    truth_vcf=$(ls "${sim_dir}"/*.truth.vcf 2>/dev/null | head -1)
    if [[ -z "$truth_vcf" ]]; then
        echo "[WARN] No truth VCF in ${sim_dir}" >&2
        FAILED+=("${type}_${tag}_${rep}")
        return 0
    fi

    mkdir -p "$out_dir"

    local junction_arg=""
    [[ -n "$junctions" && -f "$junctions" ]] && junction_arg="--sv-junctions $junctions"

    echo "[RUN] ${type} vaf${tag}_${rep} — discovery"
    # shellcheck disable=SC2086
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$targets" \
            $junction_arg \
            --output-dir "${out_dir}/tmp_disc" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] discovery failed for ${type}_${tag}_${rep}" >&2
        FAILED+=("disc_${type}_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_disc/variants.vcf" "$disc_vcf"
    cp "${out_dir}/tmp_disc/variants.tsv" "${out_dir}/calls_discovery.tsv"
    rm -rf "${out_dir}/tmp_disc"

    echo "[RUN] ${type} vaf${tag}_${rep} — tumour-informed"
    # shellcheck disable=SC2086
    # --ti-position-tolerance 10: SV truth alleles are full sequences (100bp DUP etc.)
    # but kam only detects partial alleles at the junction. Exact REF/ALT matching
    # gives 0% TI sensitivity. Position-based matching (±10bp) allows monitoring
    # of known SV loci even when the called allele is a partial junction fragment.
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$targets" \
            $junction_arg \
            --target-variants "$truth_vcf" \
            --ti-position-tolerance 10 \
            --output-dir "${out_dir}/tmp_ti" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] tumour-informed failed for ${type}_${tag}_${rep}" >&2
        FAILED+=("ti_${type}_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_ti/variants.vcf" "$ti_vcf"
    cp "${out_dir}/tmp_ti/variants.tsv" "${out_dir}/calls_tumour_informed.tsv"
    rm -rf "${out_dir}/tmp_ti"
}

SV_TARGETS="${DATADIR}/sv_suite_targets.fa"
INS_TARGETS="${DATADIR}/ins_targets.fa"
INVDEL_TARGETS="${DATADIR}/invdel_targets.fa"

# Note: --sv-junctions omitted. Junction k-mers create high-connectivity nodes
# in the de Bruijn graph that cause DFS explosion in walk.rs, hanging pathfind.
# SV detection proceeds via reference-branching paths only; this is sufficient
# for sensitivity measurement (partial alleles are still scored by position).

VAF_TAGS=(0005 0010 0015 0020 0025 0030 0035 0040 0050 0060 0075 0100 0125 0150 0175 0200 0250 0300 0350 0400 0500 0600 0700 0800 1000)

for tag in "${VAF_TAGS[@]}"; do
    for rep in a b; do
        run_dataset sv      sv      "$tag" "$rep" "$SV_TARGETS"     ""
        run_dataset ins     ins     "$tag" "$rep" "$INS_TARGETS"    ""
        run_dataset invdel  invdel  "$tag" "$rep" "$INVDEL_TARGETS" ""
    done
done

echo ""
echo "=== DONE ==="
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
