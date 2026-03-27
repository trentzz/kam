#!/usr/bin/env bash
# run.sh — Tumour-informed pipeline wrapper for kam
#
# Usage:
#   ./run.sh R1 R2 REFERENCE OUTPUT_DIR MODE
#
#   R1          Path to plasma R1 FASTQ (may be .gz)
#   R2          Path to plasma R2 FASTQ (may be .gz)
#   REFERENCE   Path to reference FASTA (hg38.fa or similar; .fai built automatically)
#   OUTPUT_DIR  Directory to write kam results
#   MODE        One of: discovery | tumour-informed | tumour-informed-altseq | alt-seq
#
# Examples:
#   ./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ discovery
#   ./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ tumour-informed
#   ./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ tumour-informed-altseq
#   ./run.sh plasma_R1.fq.gz plasma_R2.fq.gz hg38.fa results/ alt-seq
#
# Override tool paths by setting KAM or MULTISEQEX before running:
#   KAM=/usr/local/bin/kam ./run.sh ...
#   MULTISEQEX=~/.cargo/bin/multiseqex ./run.sh ...

set -euo pipefail

KAM="${KAM:-kam}"
MULTISEQEX="${MULTISEQEX:-multiseqex}"

# Resolve the directory containing this script so configs and inputs can be
# located relative to it regardless of the caller's working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIGS="${SCRIPT_DIR}/configs"
INPUTS="${SCRIPT_DIR}/inputs"

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 R1 R2 REFERENCE OUTPUT_DIR MODE"
    echo "MODE: discovery | tumour-informed | tumour-informed-altseq | alt-seq"
    exit 1
fi

R1="$1"
R2="$2"
REFERENCE="$3"
OUTPUT_DIR="$4"
MODE="$5"

mkdir -p "${OUTPUT_DIR}"

case "${MODE}" in

    discovery)
        echo "=== Discovery mode ==="

        echo "[1/2] Extracting panel targets from BED..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --bed "${INPUTS}/panel.bed" \
            --name-template "{chr}:{start}-{end}" \
            -o "${OUTPUT_DIR}/targets.fa"

        echo "[2/2] Running kam (discovery)..."
        "${KAM}" run \
            --config "${CONFIGS}/discovery.toml" \
            --r1 "${R1}" \
            --r2 "${R2}" \
            --targets "${OUTPUT_DIR}/targets.fa" \
            --output-dir "${OUTPUT_DIR}"

        echo "Done. Results in ${OUTPUT_DIR}/"
        ;;

    tumour-informed)
        echo "=== Tumour-informed mode (normal) ==="

        echo "[1/2] Extracting per-mutation target windows from VCF..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --vcf "${INPUTS}/mutations.vcf" \
            --flank 200 \
            --name-template "{chr}:{start}-{end}" \
            -o "${OUTPUT_DIR}/targets.fa"

        echo "[2/2] Running kam (tumour-informed)..."
        "${KAM}" run \
            --config "${CONFIGS}/tumour-informed.toml" \
            --r1 "${R1}" \
            --r2 "${R2}" \
            --targets "${OUTPUT_DIR}/targets.fa" \
            --target-variants "${INPUTS}/mutations.vcf" \
            --output-dir "${OUTPUT_DIR}"

        echo "Done. Results in ${OUTPUT_DIR}/"
        ;;

    tumour-informed-altseq)
        echo "=== Tumour-informed mode (alt-seq k-mer boost) ==="

        echo "[1/3] Extracting per-mutation target windows from VCF..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --vcf "${INPUTS}/mutations.vcf" \
            --flank 200 \
            --name-template "{chr}:{start}-{end}" \
            -o "${OUTPUT_DIR}/targets.fa"

        echo "[2/3] Generating alt allele sequences for k-mer allowlist..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --vcf "${INPUTS}/mutations.vcf" \
            --alt-seq \
            -o "${OUTPUT_DIR}/alt_seqs.fa"

        echo "[3/3] Running kam (tumour-informed, alt-seq boost)..."
        "${KAM}" run \
            --config "${CONFIGS}/tumour-informed-altseq.toml" \
            --r1 "${R1}" \
            --r2 "${R2}" \
            --targets "${OUTPUT_DIR}/targets.fa" \
            --sv-junctions "${OUTPUT_DIR}/alt_seqs.fa" \
            --target-variants "${INPUTS}/mutations.vcf" \
            --output-dir "${OUTPUT_DIR}"

        echo "Done. Results in ${OUTPUT_DIR}/"
        ;;

    alt-seq)
        echo "=== Alt-seq / SV junction mode ==="

        echo "[1/3] Extracting SV reference windows from VCF..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --vcf "${INPUTS}/sv_mutations.vcf" \
            --flank 200 \
            --name-template "{chr}:{start}-{end}" \
            -o "${OUTPUT_DIR}/sv_targets.fa"

        echo "[2/3] Generating SV junction sequences (alt allele contexts)..."
        "${MULTISEQEX}" "${REFERENCE}" \
            --vcf "${INPUTS}/sv_mutations.vcf" \
            --alt-seq \
            --flank 100 \
            -o "${OUTPUT_DIR}/sv_junctions.fa"

        echo "[3/3] Running kam (SV detection, tumour-informed)..."
        "${KAM}" run \
            --config "${CONFIGS}/alt-seq.toml" \
            --r1 "${R1}" \
            --r2 "${R2}" \
            --targets "${OUTPUT_DIR}/sv_targets.fa" \
            --sv-junctions "${OUTPUT_DIR}/sv_junctions.fa" \
            --target-variants "${INPUTS}/sv_mutations.vcf" \
            --output-dir "${OUTPUT_DIR}"

        echo "Done. Results in ${OUTPUT_DIR}/"
        ;;

    *)
        echo "Unknown mode: ${MODE}"
        echo "MODE must be one of: discovery | tumour-informed | tumour-informed-altseq | alt-seq"
        exit 1
        ;;
esac
