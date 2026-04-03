#!/usr/bin/env bash
# Run kam on the COLO829 somatic SV benchmark.
#
# IMPORTANT: COLO829 data has no UMI.
#   Each read pair is treated as a single molecule.
#   This requires kam's --no-umi mode, which is NOT yet implemented.
#   When that mode becomes available, remove the exit below and run
#   the full pipeline. Until then, this script documents the intended
#   invocation and performs preflight checks only.
#
# Dataset:
#   Tumour:  ERR1341793 — COLO829 Illumina HiSeq X Ten WGS (~100 GB)
#   Normal:  ERR1341794 — COLO829BL Illumina WGS (~80 GB)
#   Truth:   COLO829_truthset_somatic_v4.1.vcf (GRCh37, 68 validated SVs)
#            COLO829_truthset_somatic_hg38.vcf  (hg38 liftover)
#
# Truth set:
#   68 somatic SVs: 38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS.
#   Confirmed by five technologies (Illumina, ONT, PacBio, 10x, Bionano).
#   Source: Zenodo record 4716169, van Dijk et al., Cell Genomics 2022.
#
# Evaluation:
#   Use Truvari for SV comparison (see seqc2/README.md for rtg vcfeval notes).
#   Example:
#     truvari bench \
#         --base data/COLO829_truthset_somatic_v4.1.vcf \
#         --comp results/kam_colo829_sv/calls_discovery.vcf.gz \
#         --output results/truvari_discovery/ \
#         --refdist 1000 --pctsim 0.7 --pctsize 0.7
#
# Prerequisites:
#   - Download data:         sh docs/benchmarking/05-public/scripts/download_colo829.sh
#   - Design SV targets:     python3 design_sv_targets.py (in this directory)
#   - Build kam:             cargo build --release
#
# Usage:
#   bash run_kam.sh [--force]
#
# --force: re-run even if output already exists.
#
# Outputs are written to:
#   docs/benchmarking/05-public/colo829/results/

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
KAM="${REPO}/target/release/kam"
SCRIPTDIR="${REPO}/docs/benchmarking/05-public/colo829"
DATADIR="${SCRIPTDIR}/data"
RESDIR="${SCRIPTDIR}/results"
CONFIG="${RESDIR}/config.toml"
FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

# ── No-UMI mode gate ──────────────────────────────────────────────────────────

# Remove this block when --no-umi is implemented in kam.
echo "[BLOCKED] COLO829 requires --no-umi mode, which is not yet implemented." >&2
echo "          This script documents the intended run configuration." >&2
echo "          Track progress in docs/benchmarking/05-public/colo829/README.md." >&2
echo "" >&2
echo "          When --no-umi is available, re-run this script." >&2
exit 1

# ── Preflight checks ──────────────────────────────────────────────────────────
# (Unreachable until the gate above is removed.)

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    echo "        Run: cargo build --release" >&2
    exit 1
fi

TUMOUR_R1="${DATADIR}/ERR1341793_1.fastq.gz"
TUMOUR_R2="${DATADIR}/ERR1341793_2.fastq.gz"
TRUTH_VCF="${DATADIR}/COLO829_truthset_somatic_v4.1.vcf"
TARGETS="${DATADIR}/colo829_sv_targets.fa"
JUNCTIONS="${DATADIR}/colo829_sv_junctions.fa"

for f in "$TUMOUR_R1" "$TUMOUR_R2"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] FASTQ not found: $f" >&2
        echo "        Run: sh ${REPO}/docs/benchmarking/05-public/scripts/download_colo829.sh" >&2
        exit 1
    fi
done

if [[ ! -f "$TRUTH_VCF" ]]; then
    echo "[ERROR] Truth VCF not found: $TRUTH_VCF" >&2
    echo "        Run: sh ${REPO}/docs/benchmarking/05-public/scripts/download_colo829.sh" >&2
    exit 1
fi

if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] SV target sequences not found: $TARGETS" >&2
    echo "        Generate with:" >&2
    echo "          python3 ${SCRIPTDIR}/scripts/design_sv_targets.py \\" >&2
    echo "              --vcf $TRUTH_VCF \\" >&2
    echo "              --ref <path/to/GRCh37.fa> \\" >&2
    echo "              --targets-out $TARGETS \\" >&2
    echo "              --junctions-out $JUNCTIONS" >&2
    exit 1
fi

OUT="${RESDIR}/kam_colo829_sv"
DONE_DISC="${OUT}/calls_discovery.vcf"
DONE_TI="${OUT}/calls_tumour_informed.vcf"

if [[ "$FORCE" == "false" && -f "$DONE_DISC" && -f "$DONE_TI" ]]; then
    echo "[SKIP] Outputs already exist. Use --force to re-run."
    exit 0
fi

# ── Write config ──────────────────────────────────────────────────────────────

mkdir -p "$RESDIR" "$OUT"

cat > "$CONFIG" << 'TOML'
# kam config for COLO829 somatic SV benchmark.
#
# COLO829 has no UMI. umi_length = 0 and skip_length = 0.
# Each read pair is treated as one molecule (no UMI grouping).
# duplex = false: no duplex information is available.
#
# SV calling settings:
#   sv_min_confidence = 0.95 — require high posterior for SV PASS calls.
#   sv_min_alt_molecules = 3 — require at least 3 supporting molecules.
#   sv_strand_bias_threshold = 1.0 — disabled; inversion reads are
#     structurally strand-biased and the SNV threshold is inappropriate.
#   ti_position_tolerance = 500 — SV truth alleles are full genomic
#     sequences; kam detects partial junction alleles. Position-based
#     matching allows monitoring of known SV loci.

[chemistry]
umi_length  = 0
skip_length = 0
duplex      = false

[assembly]
min_family_size = 1

[indexing]
kmer_size = 31

[calling]
min_confidence        = 0.95
min_alt_molecules     = 3
sv_min_confidence     = 0.95
sv_min_alt_molecules  = 3
sv_strand_bias_threshold = 1.0
max_vaf               = 0.80

[output]
output_format = "tsv,vcf"
TOML

echo "[INFO] Config written: $CONFIG"

# ── Run kam: discovery mode ───────────────────────────────────────────────────

# Determine whether junction sequences are available.
JUNCTION_ARG=""
if [[ -f "$JUNCTIONS" ]]; then
    JUNCTION_ARG="--sv-junctions $JUNCTIONS"
    echo "[INFO] SV junction sequences: $JUNCTIONS"
else
    echo "[INFO] No junction FASTA found; running without --sv-junctions."
fi

echo "[RUN] kam run — discovery mode"
# shellcheck disable=SC2086
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS" \
    $JUNCTION_ARG \
    --output-dir "${OUT}/tmp_disc" \
    --qc-output "${OUT}/qc_disc.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT}/tmp_disc/variants.vcf" "${OUT}/calls_discovery.vcf"
cp "${OUT}/tmp_disc/variants.tsv" "${OUT}/calls_discovery.tsv"
rm -rf "${OUT}/tmp_disc"

# ── Run kam: tumour-informed mode ─────────────────────────────────────────────

echo "[RUN] kam run — tumour-informed mode"
# --ti-position-tolerance 500: SV truth alleles are full-length genomic
# sequences that kam cannot exactly reproduce as k-mer alleles.
# Position-based matching within ±500 bp allows tumour-informed monitoring
# of known SV breakpoints.
# shellcheck disable=SC2086
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS" \
    $JUNCTION_ARG \
    --target-variants "$TRUTH_VCF" \
    --ti-position-tolerance 500 \
    --output-dir "${OUT}/tmp_ti" \
    --qc-output "${OUT}/qc_ti.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT}/tmp_ti/variants.vcf" "${OUT}/calls_tumour_informed.vcf"
cp "${OUT}/tmp_ti/variants.tsv" "${OUT}/calls_tumour_informed.tsv"
rm -rf "${OUT}/tmp_ti"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "=== DONE ==="
echo ""
echo "Outputs:"
ls -lh "${OUT}/" 2>/dev/null || echo "  (none)"
echo ""
echo "Truth set: 68 validated somatic SVs (38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS)"
echo ""
echo "Evaluate with Truvari:"
echo "  bgzip ${OUT}/calls_discovery.vcf"
echo "  tabix ${OUT}/calls_discovery.vcf.gz"
echo "  truvari bench \\"
echo "      --base ${TRUTH_VCF} \\"
echo "      --comp ${OUT}/calls_discovery.vcf.gz \\"
echo "      --output ${RESDIR}/truvari_discovery/ \\"
echo "      --refdist 1000 --pctsim 0.7 --pctsize 0.7"
echo ""
echo "See docs/benchmarking/05-public/colo829/README.md for full evaluation notes."
