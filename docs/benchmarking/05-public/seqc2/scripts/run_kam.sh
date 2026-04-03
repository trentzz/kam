#!/usr/bin/env bash
# Run kam on the SEQC2 HCC1395 somatic mutation benchmark.
#
# IMPORTANT: SEQC2 data has no UMI.
#   Each read pair is treated as a single molecule.
#   This requires kam's --no-umi mode, which is NOT yet implemented.
#   When that mode becomes available, remove the exit below and run
#   the full pipeline. Until then, this script documents the intended
#   invocation and performs preflight checks only.
#
# Dataset:
#   Tumour:  SRR7890824 — HCC1395 Illumina WGS (~50 GB)
#   Normal:  SRR7890827 — HCC1395BL Illumina WGS (~50 GB)
#   Truth:   HCC1395_truth_SNV.vcf.gz, HCC1395_truth_indel.vcf.gz
#
# Truth set evaluation:
#   Use rtg vcfeval or hap.py. Restrict to PASS variants only.
#   See docs/benchmarking/05-public/seqc2/README.md.
#
# Prerequisites:
#   - Download data:   sh docs/benchmarking/05-public/scripts/download_seqc2.sh
#   - Design targets:  python3 design_targets.py (in this directory)
#   - Build kam:       cargo build --release
#
# Usage:
#   bash run_kam.sh [--force]
#
# --force: re-run even if output already exists.
#
# Outputs are written to:
#   docs/benchmarking/05-public/seqc2/results/

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
KAM="${REPO}/target/release/kam"
SCRIPTDIR="${REPO}/docs/benchmarking/05-public/seqc2"
DATADIR="${SCRIPTDIR}/data"
RESDIR="${SCRIPTDIR}/results"
CONFIG="${RESDIR}/config.toml"
FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

# ── No-UMI mode gate ──────────────────────────────────────────────────────────

# Remove this block when --no-umi is implemented in kam.
echo "[BLOCKED] SEQC2 requires --no-umi mode, which is not yet implemented." >&2
echo "          This script documents the intended run configuration." >&2
echo "          Track progress in docs/benchmarking/05-public/seqc2/README.md." >&2
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

TUMOUR_R1="${DATADIR}/SRR7890824_1.fastq.gz"
TUMOUR_R2="${DATADIR}/SRR7890824_2.fastq.gz"
NORMAL_R1="${DATADIR}/SRR7890827_1.fastq.gz"
NORMAL_R2="${DATADIR}/SRR7890827_2.fastq.gz"
TRUTH_SNV="${DATADIR}/HCC1395_truth_SNV.vcf.gz"
TRUTH_INDEL="${DATADIR}/HCC1395_truth_indel.vcf.gz"
TARGETS="${DATADIR}/seqc2_targets.fa"

for f in "$TUMOUR_R1" "$TUMOUR_R2" "$NORMAL_R1" "$NORMAL_R2"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] FASTQ not found: $f" >&2
        echo "        Run: sh ${REPO}/docs/benchmarking/05-public/scripts/download_seqc2.sh" >&2
        exit 1
    fi
done

for f in "$TRUTH_SNV" "$TRUTH_INDEL"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] Truth VCF not found: $f" >&2
        echo "        Run: sh ${REPO}/docs/benchmarking/05-public/scripts/download_seqc2.sh" >&2
        exit 1
    fi
done

if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Target sequences FASTA not found: $TARGETS" >&2
    echo "        Generate it with:" >&2
    echo "          python3 ${SCRIPTDIR}/scripts/design_targets.py \\" >&2
    echo "              --vcf $TRUTH_SNV \\" >&2
    echo "              --ref <path/to/GRCh38.fa> \\" >&2
    echo "              --output $TARGETS" >&2
    exit 1
fi

OUT_SNV="${RESDIR}/kam_seqc2_snv"
OUT_INDEL="${RESDIR}/kam_seqc2_indel"
DONE_SNV="${OUT_SNV}/calls_discovery.vcf"
DONE_INDEL="${OUT_INDEL}/calls_discovery.vcf"

if [[ "$FORCE" == "false" && -f "$DONE_SNV" && -f "$DONE_INDEL" ]]; then
    echo "[SKIP] Outputs already exist. Use --force to re-run."
    exit 0
fi

# ── Write config ──────────────────────────────────────────────────────────────

mkdir -p "$RESDIR" "$OUT_SNV" "$OUT_INDEL"

cat > "$CONFIG" << 'TOML'
# kam config for SEQC2 HCC1395 benchmark.
#
# SEQC2 has no UMI. umi_length = 0 and skip_length = 0.
# Each read pair is treated as one molecule (no UMI grouping).
# duplex = false: no duplex information is available.
#
# This run benchmarks the variant calling component of kam only.
# Molecule assembly is trivial without UMI deduplication.

[chemistry]
umi_length  = 0
skip_length = 0
duplex      = false

[assembly]
min_family_size = 1

[indexing]
kmer_size = 31

[calling]
min_confidence    = 0.95
min_alt_molecules = 3
# max_vaf suppresses germline heterozygous variants (VAF ~0.5).
max_vaf = 0.35

[output]
output_format = "tsv,vcf"
TOML

echo "[INFO] Config written: $CONFIG"

# ── Run kam: SNV targets ──────────────────────────────────────────────────────

TARGETS_SNV="${DATADIR}/seqc2_snv_targets.fa"
TARGETS_INDEL="${DATADIR}/seqc2_indel_targets.fa"

echo "[RUN] kam run — SNV targets, discovery mode"
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS_SNV" \
    --output-dir "${OUT_SNV}/tmp" \
    --qc-output "${OUT_SNV}/qc.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT_SNV}/tmp/variants.vcf" "${OUT_SNV}/calls_discovery.vcf"
cp "${OUT_SNV}/tmp/variants.tsv" "${OUT_SNV}/calls_discovery.tsv"
rm -rf "${OUT_SNV}/tmp"

echo "[RUN] kam run — SNV targets, tumour-informed mode"
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS_SNV" \
    --target-variants "$TRUTH_SNV" \
    --output-dir "${OUT_SNV}/tmp_ti" \
    --qc-output "${OUT_SNV}/qc_ti.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT_SNV}/tmp_ti/variants.vcf" "${OUT_SNV}/calls_tumour_informed.vcf"
cp "${OUT_SNV}/tmp_ti/variants.tsv" "${OUT_SNV}/calls_tumour_informed.tsv"
rm -rf "${OUT_SNV}/tmp_ti"

# ── Run kam: indel targets ────────────────────────────────────────────────────

echo "[RUN] kam run — indel targets, discovery mode"
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS_INDEL" \
    --output-dir "${OUT_INDEL}/tmp" \
    --qc-output "${OUT_INDEL}/qc.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT_INDEL}/tmp/variants.vcf" "${OUT_INDEL}/calls_discovery.vcf"
cp "${OUT_INDEL}/tmp/variants.tsv" "${OUT_INDEL}/calls_discovery.tsv"
rm -rf "${OUT_INDEL}/tmp"

echo "[RUN] kam run — indel targets, tumour-informed mode"
"$KAM" run \
    --config "$CONFIG" \
    --r1 "$TUMOUR_R1" \
    --r2 "$TUMOUR_R2" \
    --targets "$TARGETS_INDEL" \
    --target-variants "$TRUTH_INDEL" \
    --output-dir "${OUT_INDEL}/tmp_ti" \
    --qc-output "${OUT_INDEL}/qc_ti.json" \
    2>&1 | grep -E "^\[|ERROR" || true

cp "${OUT_INDEL}/tmp_ti/variants.vcf" "${OUT_INDEL}/calls_tumour_informed.vcf"
cp "${OUT_INDEL}/tmp_ti/variants.tsv" "${OUT_INDEL}/calls_tumour_informed.tsv"
rm -rf "${OUT_INDEL}/tmp_ti"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "=== DONE ==="
echo ""
echo "Outputs:"
ls -lh "${OUT_SNV}/" "${OUT_INDEL}/" 2>/dev/null || echo "  (none)"
echo ""
echo "Evaluate against truth VCF:"
echo "  rtg vcfeval \\"
echo "      --baseline ${TRUTH_SNV} \\"
echo "      --calls ${OUT_SNV}/calls_discovery.vcf.gz \\"
echo "      --template GRCh38.sdf \\"
echo "      --output ${RESDIR}/rtg_eval_snv_discovery/"
echo ""
echo "Restrict evaluation to PASS variants in the truth VCF."
echo "See docs/benchmarking/05-public/seqc2/README.md for full evaluation steps."
