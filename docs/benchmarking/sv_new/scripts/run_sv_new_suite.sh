#!/usr/bin/env bash
# Run kam on all sv_new benchmark datasets: InvDel, NovIns, and Fusion.
#
# For each varforge-simulated dataset:
#   1. Discovery mode  — no --target-variants
#   2. Tumour-informed mode — --target-variants using the varforge truth VCF
#
# Fusion datasets require a mixing step: two separate varforge runs (fusion
# reads + wild-type reads) are concatenated before kam is invoked.
#
# Outputs per dataset:
#   results/kam_{type}_vaf{tag}_{rep}/calls_discovery.vcf
#   results/kam_{type}_vaf{tag}_{rep}/calls_discovery.tsv
#   results/kam_{type}_vaf{tag}_{rep}/calls_tumour_informed.vcf
#   results/kam_{type}_vaf{tag}_{rep}/calls_tumour_informed.tsv
#
# Usage:
#   bash run_sv_new_suite.sh [--force]
#
#   --force  Re-run even if output files already exist.
#
# Requires: target/release/kam, varforge, Python 3.
# Run from the repo root, or from any directory (paths are resolved from the
# script location).

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../" && pwd)"
KAM="${REPO}/target/release/kam"
EXPDIR="${REPO}/docs/benchmarking/sv_new"
CFGDIR="${EXPDIR}/configs"
DATADIR="${EXPDIR}/data"
RESDIR="${EXPDIR}/results"
SV_DATADIR="${REPO}/docs/benchmarking/sv/data"

INVDEL_TARGETS="${SV_DATADIR}/invdel_targets.fa"
NOVINS_TARGETS="${SV_DATADIR}/ins_targets.fa"
# Fusion uses the standard ref as primary targets and adds the fusion target
# via --fusion-targets.
FUSION_TARGETS="${SV_DATADIR}/ref.fa"
FUSION_JUNCTION_TARGETS="${DATADIR}/fusion_targets.fa"

FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    exit 1
fi

if ! command -v varforge >/dev/null 2>&1; then
    echo "[ERROR] varforge not found in PATH" >&2
    exit 1
fi

FAILED=()

# ── Non-fusion datasets (invdel, novins) ──────────────────────────────────────

run_standard_dataset() {
    # Run varforge simulation and both kam modes for a standard (non-fusion) type.
    #
    # Arguments:
    #   type    — invdel or novins
    #   tag     — zero-padded 4-digit VAF tag, e.g. 0100
    #   rep     — replicate letter: a or b
    #   targets — path to the kam targets FASTA
    local type="$1"
    local tag="$2"
    local rep="$3"
    local targets="$4"

    local cfg="${CFGDIR}/${type}_vaf${tag}_${rep}.yaml"
    local sim_dir="${RESDIR}/sim_${type}_vaf${tag}_${rep}"
    local out_dir="${RESDIR}/kam_${type}_vaf${tag}_${rep}"

    local disc_vcf="${out_dir}/calls_discovery.vcf"
    local ti_vcf="${out_dir}/calls_tumour_informed.vcf"

    if [[ "$FORCE" == "false" && -f "$disc_vcf" && -f "$ti_vcf" ]]; then
        return 0
    fi

    # Run varforge if the simulation output does not already exist.
    if [[ "$FORCE" == "true" || ! -d "$sim_dir" ]]; then
        echo "[SIM] ${type} vaf${tag}_${rep}"
        varforge simulate -c "$cfg"
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

    echo "[RUN] ${type} vaf${tag}_${rep} — discovery"
    # shellcheck disable=SC2086
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$targets" \
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
    # --ti-position-tolerance 10: SV alleles reported by kam are partial junction
    # fragments. Exact REF/ALT matching gives 0% tumour-informed sensitivity for
    # large SVs. Position-based matching (±10 bp) allows monitoring known loci.
    # shellcheck disable=SC2086
    if ! "$KAM" run \
            --r1 "$r1" --r2 "$r2" \
            --targets "$targets" \
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

# ── Fusion datasets ───────────────────────────────────────────────────────────

run_fusion_dataset() {
    # Fusions require two varforge runs (fusion + wild-type) followed by
    # FASTQ concatenation to produce the mixed sample.
    #
    # Arguments:
    #   tag  — zero-padded 4-digit VAF tag, e.g. 0100
    #   rep  — replicate letter: a or b
    local tag="$1"
    local rep="$2"

    local fusion_cfg="${CFGDIR}/fusion_vaf${tag}_${rep}.yaml"
    local wt_cfg="${CFGDIR}/fusion_vaf${tag}_${rep}_wt.yaml"

    local fusion_sim="${RESDIR}/sim_fusion_vaf${tag}_${rep}"
    local wt_sim="${RESDIR}/sim_fusion_vaf${tag}_${rep}_wt"
    local mixed_dir="${RESDIR}/sim_fusion_vaf${tag}_${rep}_mixed"

    local out_dir="${RESDIR}/kam_fusion_vaf${tag}_${rep}"
    local disc_vcf="${out_dir}/calls_discovery.vcf"
    local ti_vcf="${out_dir}/calls_tumour_informed.vcf"

    if [[ "$FORCE" == "false" && -f "$disc_vcf" && -f "$ti_vcf" ]]; then
        return 0
    fi

    # Step 1: simulate fusion reads.
    if [[ "$FORCE" == "true" || ! -d "$fusion_sim" ]]; then
        echo "[SIM] fusion vaf${tag}_${rep} — fusion reads"
        varforge simulate -c "$fusion_cfg"
    fi

    # Step 2: simulate wild-type background reads.
    if [[ "$FORCE" == "true" || ! -d "$wt_sim" ]]; then
        echo "[SIM] fusion vaf${tag}_${rep} — wild-type reads"
        varforge simulate -c "$wt_cfg"
    fi

    # Step 3: mix by concatenating the two FASTQ sets.
    mkdir -p "$mixed_dir"

    local fusion_r1 fusion_r2 wt_r1 wt_r2
    fusion_r1=$(ls "${fusion_sim}"/*_R1.fastq.gz 2>/dev/null | head -1)
    fusion_r2=$(ls "${fusion_sim}"/*_R2.fastq.gz 2>/dev/null | head -1)
    wt_r1=$(ls "${wt_sim}"/*_R1.fastq.gz 2>/dev/null | head -1)
    wt_r2=$(ls "${wt_sim}"/*_R2.fastq.gz 2>/dev/null | head -1)

    if [[ -z "$fusion_r1" || -z "$fusion_r2" || -z "$wt_r1" || -z "$wt_r2" ]]; then
        echo "[WARN] Missing FASTQs for fusion vaf${tag}_${rep}" >&2
        FAILED+=("fusion_${tag}_${rep}")
        return 0
    fi

    local mixed_r1="${mixed_dir}/sample_R1.fastq.gz"
    local mixed_r2="${mixed_dir}/sample_R2.fastq.gz"

    if [[ "$FORCE" == "true" || ! -f "$mixed_r1" ]]; then
        echo "[MIX] fusion vaf${tag}_${rep} — concatenating FASTQs"
        cat "$fusion_r1" "$wt_r1" > "$mixed_r1"
        cat "$fusion_r2" "$wt_r2" > "$mixed_r2"
    fi

    # Locate the truth VCF from the sv_new data directory.
    local truth_vcf="${DATADIR}/truth_fusion_vaf${tag}_${rep}.vcf"
    if [[ ! -f "$truth_vcf" ]]; then
        echo "[WARN] No truth VCF: ${truth_vcf}" >&2
        FAILED+=("fusion_${tag}_${rep}")
        return 0
    fi

    mkdir -p "$out_dir"

    echo "[RUN] fusion vaf${tag}_${rep} — discovery"
    if ! "$KAM" run \
            --r1 "$mixed_r1" --r2 "$mixed_r2" \
            --targets "$FUSION_TARGETS" \
            --fusion-targets "$FUSION_JUNCTION_TARGETS" \
            --output-dir "${out_dir}/tmp_disc" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] discovery failed for fusion_${tag}_${rep}" >&2
        FAILED+=("disc_fusion_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_disc/variants.vcf" "$disc_vcf"
    cp "${out_dir}/tmp_disc/variants.tsv" "${out_dir}/calls_discovery.tsv"
    rm -rf "${out_dir}/tmp_disc"

    echo "[RUN] fusion vaf${tag}_${rep} — tumour-informed"
    if ! "$KAM" run \
            --r1 "$mixed_r1" --r2 "$mixed_r2" \
            --targets "$FUSION_TARGETS" \
            --fusion-targets "$FUSION_JUNCTION_TARGETS" \
            --target-variants "$truth_vcf" \
            --ti-position-tolerance 10 \
            --output-dir "${out_dir}/tmp_ti" \
            --output-format vcf,tsv \
            2>&1 | grep -E "^\[run\]|ERROR"; then
        echo "[ERROR] tumour-informed failed for fusion_${tag}_${rep}" >&2
        FAILED+=("ti_fusion_${tag}_${rep}")
        return 0
    fi
    cp "${out_dir}/tmp_ti/variants.vcf" "$ti_vcf"
    cp "${out_dir}/tmp_ti/variants.tsv" "${out_dir}/calls_tumour_informed.tsv"
    rm -rf "${out_dir}/tmp_ti"
}

# ── Main loop ─────────────────────────────────────────────────────────────────

VAF_TAGS=(0005 0010 0015 0020 0025 0030 0035 0040 0050 0060 0075 0100 0125 0150 0175 0200 0250 0300 0350 0400 0500 0600 0700 0800 1000)

for tag in "${VAF_TAGS[@]}"; do
    for rep in a b; do
        run_standard_dataset invdel  "$tag" "$rep" "$INVDEL_TARGETS"
        run_standard_dataset novins  "$tag" "$rep" "$NOVINS_TARGETS"
        run_fusion_dataset           "$tag" "$rep"
    done
done

echo ""
echo "[BUILD] Rebuilding per-sample directories..."
python3 "${REPO}/docs/benchmarking/build_sample_dirs.py"

echo ""
echo "=== DONE ==="
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
