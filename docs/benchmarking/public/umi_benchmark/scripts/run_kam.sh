#!/usr/bin/env bash
# Run kam on the UMI Clustering Benchmark dataset (SRR6794144).
#
# This script validates kam's molecule assembly on real cfDNA UMI data.
# There is no SNV/indel truth set for this accession. The run is
# discovery-only: the goal is to confirm that kam groups reads by UMI
# correctly, not to measure variant calling accuracy.
#
# Dataset:
#   SRR6794144 — cfDNA panel sequencing, 12 bp simplex UMI in read name.
#   Chemistry: NOT Twist 5M2S+T. The UMI is embedded in the read name
#   and must be extracted before running kam (see README.md).
#
# Prerequisites:
#   - Download data with: sh docs/benchmarking/public/scripts/download_umi_benchmark.sh
#   - Extract the UMI from the read name with umi_tools extract (see below).
#   - Build kam: cargo build --release
#
# UMI extraction (required before running this script):
#   umi_tools extract \
#       --bc-pattern=NNNNNNNNNNNN \
#       --stdin data/SRR6794144_1.fastq.gz \
#       --stdout data/SRR6794144_1.umi.fastq.gz \
#       --read2-in data/SRR6794144_2.fastq.gz \
#       --read2-out data/SRR6794144_2.umi.fastq.gz \
#       --log data/extract.log
#
# After extraction, the UMI occupies the first 12 bases of R1 with no
# spacer. The config below sets umi_length=12 and skip_length=0.
#
# Usage:
#   bash run_kam.sh [--force]
#
# --force: re-run even if output already exists.
#
# Outputs are written to:
#   docs/benchmarking/public/umi_benchmark/results/

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
KAM="${REPO}/target/release/kam"
SCRIPTDIR="${REPO}/docs/benchmarking/public/umi_benchmark"
DATADIR="${SCRIPTDIR}/data"
RESDIR="${SCRIPTDIR}/results"
CONFIG="${RESDIR}/config.toml"
FORCE=false

for arg in "$@"; do
    [[ "$arg" == "--force" ]] && FORCE=true
done

# ── Preflight checks ──────────────────────────────────────────────────────────

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    echo "        Run: cargo build --release" >&2
    exit 1
fi

R1="${DATADIR}/SRR6794144_1.umi.fastq.gz"
R2="${DATADIR}/SRR6794144_2.umi.fastq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[ERROR] UMI-extracted FASTQ not found." >&2
    echo "        Expected:" >&2
    echo "          $R1" >&2
    echo "          $R2" >&2
    echo "" >&2
    echo "        Step 1: Download the raw FASTQ:" >&2
    echo "          sh ${REPO}/docs/benchmarking/public/scripts/download_umi_benchmark.sh" >&2
    echo "" >&2
    echo "        Step 2: Extract the 12 bp UMI from the read name:" >&2
    echo "          umi_tools extract \\" >&2
    echo "              --bc-pattern=NNNNNNNNNNNN \\" >&2
    echo "              --stdin ${DATADIR}/SRR6794144_1.fastq.gz \\" >&2
    echo "              --stdout ${R1} \\" >&2
    echo "              --read2-in ${DATADIR}/SRR6794144_2.fastq.gz \\" >&2
    echo "              --read2-out ${R2} \\" >&2
    echo "              --log ${DATADIR}/extract.log" >&2
    exit 1
fi

TARGETS="${DATADIR}/targets.fa"
if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Target sequences FASTA not found: $TARGETS" >&2
    echo "        This dataset has no associated panel design." >&2
    echo "        Provide a panel FASTA or run in de-novo mode (not yet implemented)." >&2
    echo "        For assembly-only validation, a dummy single-base target is sufficient:" >&2
    echo "          printf '>dummy\\nACGT\\n' > $TARGETS" >&2
    exit 1
fi

OUT="${RESDIR}/kam_umi_benchmark"
DONE_MARKER="${OUT}/calls_discovery.vcf"

if [[ "$FORCE" == "false" && -f "$DONE_MARKER" ]]; then
    echo "[SKIP] Output already exists: $DONE_MARKER"
    echo "       Use --force to re-run."
    exit 0
fi

# ── Write config ──────────────────────────────────────────────────────────────

mkdir -p "$RESDIR" "$OUT"

cat > "$CONFIG" << 'TOML'
# kam config for UMI Clustering Benchmark (SRR6794144).
#
# Chemistry: 12 bp simplex UMI, no spacer.
# The UMI has been moved into the first 12 bases of R1 by umi_tools extract.
# skip_length = 0 because the original data has no 2 bp spacer.
# duplex = false: this is a simplex dataset (one-strand UMI, no duplex pairing).
#
# There is no SNV/indel truth set for this accession. The purpose of
# this run is molecule assembly validation, not variant calling accuracy.

[chemistry]
umi_length  = 12
skip_length = 0
duplex      = false
min_umi_quality = 20

[assembly]
min_family_size = 2

[indexing]
kmer_size = 31

[calling]
min_confidence    = 0.99
min_alt_molecules = 2

[output]
output_format = "tsv,vcf"
TOML

echo "[INFO] Config written: $CONFIG"

# ── Run kam ───────────────────────────────────────────────────────────────────

echo "[RUN] kam run — discovery mode (no truth set)"

"$KAM" run \
    --config "$CONFIG" \
    --r1 "$R1" \
    --r2 "$R2" \
    --targets "$TARGETS" \
    --output-dir "${OUT}/tmp" \
    --qc-output "${OUT}/qc.json" \
    2>&1 | grep -E "^\[|ERROR" || true

# Copy outputs to canonical names.
if [[ -f "${OUT}/tmp/variants.vcf" ]]; then
    cp "${OUT}/tmp/variants.vcf" "${OUT}/calls_discovery.vcf"
fi
if [[ -f "${OUT}/tmp/variants.tsv" ]]; then
    cp "${OUT}/tmp/variants.tsv" "${OUT}/calls_discovery.tsv"
fi
rm -rf "${OUT}/tmp"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "=== DONE ==="
echo ""
echo "Outputs:"
ls -lh "${OUT}/" 2>/dev/null || echo "  (none)"
echo ""
echo "Notes:"
echo "  - No truth set: this run validates molecule assembly only."
echo "  - Check qc.json for molecule counts, UMI family size distribution,"
echo "    and duplicate rates."
echo "  - Variant calls in calls_discovery.vcf/tsv are exploratory only."
echo "  - See docs/benchmarking/public/umi_benchmark/README.md for context."
