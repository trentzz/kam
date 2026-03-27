#!/usr/bin/env bash
# Tumour-informed kam pipeline using multiseqex for target generation.
#
# OVERVIEW
# --------
# This script turns a tumour mutation VCF and a pair of plasma FASTQ files
# into a VCF of tracked ctDNA mutations in two steps:
#
#   1. multiseqex: extract 401 bp windows (±200 bp) around each variant in
#      the mutation VCF from the reference FASTA → targets.fa
#
#   2. kam: index those windows, walk de Bruijn paths, call variants, and
#      apply tumour-informed filtering so only the known mutations appear as PASS.
#
# REQUIREMENTS
# ------------
#   - multiseqex >= 0.2.1 (needs --vcf support; build from source if needed)
#       cargo install --path /path/to/multiseqex
#     or:
#       cargo install --git https://github.com/trentzz/multiseqex
#
#   - kam >= 0.1.0 installed in PATH, or set KAM below to the binary path.
#
# INPUT FILE FORMATS
# ------------------
# FASTQ (plasma sample, gzip or plain):
#   plasma_R1.fq.gz   plasma_R2.fq.gz
#
# Reference FASTA (must be indexed; multiseqex builds .fai automatically):
#   reference.fa
#
# Tumour mutation VCF — the set of variants to track.
# Produced by variant calling on a matched tumour biopsy or a ctDNA baseline
# run. Standard 5-column VCF is sufficient; extra columns (QUAL, FILTER, INFO,
# FORMAT, SAMPLE) are ignored by multiseqex and kam.
#
#   ##fileformat=VCFv4.2
#   #CHROM  POS       ID  REF  ALT
#   chr7    55191823  .   T    G
#   chr17   7674220   .   C    T
#   chr12   25398284  .   C    A
#
# The ID column must be "." (dot) so that multiseqex produces plain
# "chr:start-end" FASTA headers. Named IDs would prepend the name to the
# header and break kam's coordinate parser.
#
# USAGE
# -----
#   bash examples/ti-pipeline.sh \
#       plasma_R1.fq.gz plasma_R2.fq.gz \
#       tumour_mutations.vcf \
#       reference.fa \
#       results/
#
# Or set the variables below and run directly.

set -euo pipefail

# ── Inputs ─────────────────────────────────────────────────────────────────────
R1="${1:-plasma_R1.fq.gz}"
R2="${2:-plasma_R2.fq.gz}"
MUTATIONS_VCF="${3:-tumour_mutations.vcf}"
REFERENCE="${4:-reference.fa}"
OUT_DIR="${5:-results}"

# ── Tool paths ─────────────────────────────────────────────────────────────────
KAM="${KAM:-kam}"
MULTISEQEX="${MULTISEQEX:-multiseqex}"

# ── Parameters ─────────────────────────────────────────────────────────────────
# Window size around each variant position (bp each side).
FLANK=200

# ── Derived paths ──────────────────────────────────────────────────────────────
TARGETS_FA="${OUT_DIR}/targets.fa"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/ti-pipeline.toml"

# ── Setup ──────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"

# ── Step 1: Generate target windows ───────────────────────────────────────────
# multiseqex reads the reference FASTA and extracts one 401 bp window per
# variant in the mutation VCF. Output FASTA headers are in "chr:start-end"
# format (--name-template enforces this regardless of VCF ID column content).
#
# The reference FASTA does not need a pre-built .fai; multiseqex builds one
# automatically and stores it alongside the FASTA file.
echo "[step 1/2] Extracting target windows with multiseqex"
"${MULTISEQEX}" "${REFERENCE}" \
    --vcf "${MUTATIONS_VCF}" \
    --flank "${FLANK}" \
    --name-template "{chr}:{start}-{end}" \
    -o "${TARGETS_FA}"

N_TARGETS=$(grep -c '^>' "${TARGETS_FA}")
echo "  Wrote ${N_TARGETS} target sequence(s) to ${TARGETS_FA}"

# ── Step 2: Run kam in tumour-informed mode ────────────────────────────────────
# CLI flags override the TOML values, so r1, r2, targets, target_variants, and
# output_dir are set here rather than hard-coded in the config file. The config
# file carries the chemistry, assembly, calling, and tolerance settings.
echo "[step 2/2] Running kam (tumour-informed)"
"${KAM}" run \
    --config "${CONFIG}" \
    --r1 "${R1}" \
    --r2 "${R2}" \
    --targets "${TARGETS_FA}" \
    --target-variants "${MUTATIONS_VCF}" \
    --output-dir "${OUT_DIR}"

echo ""
echo "Done. PASS calls:"
grep -c '^[^#]' "${OUT_DIR}/variants.vcf" 2>/dev/null | xargs -I{} echo "  {} total records"
grep -v '^#' "${OUT_DIR}/variants.vcf" 2>/dev/null | awk -F'\t' '$7=="PASS"' | wc -l | xargs -I{} echo "  {} PASS"
echo "  Full results: ${OUT_DIR}/"
