#!/bin/sh
# Download the UMI Clustering Benchmark template dataset (SRR6794144).
#
# Source: NCBI SRA, open access.
# Paper: Benchmarking UMI clustering tools for accurate detection of
#        low-frequency variants from deep sequencing (Scientific Reports, 2025).
#        https://www.nature.com/articles/s41598-025-33128-x
#
# Expected disk space: ~20 GB uncompressed FASTQ pair (~10 GB gzip compressed).
# Runtime: 20–60 minutes depending on network speed.
#
# Usage:
#   sh download_umi_benchmark.sh [OUTPUT_DIR]
#
# If OUTPUT_DIR is not given, files are placed in ./data/umi_benchmark/.

set -eu

ACCESSION="SRR6794144"
OUTPUT_DIR="${1:-data/umi_benchmark}"

# Check for required tools.
if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "Error: fasterq-dump not found." >&2
    echo "Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit" >&2
    exit 1
fi

if ! command -v pigz >/dev/null 2>&1 && ! command -v gzip >/dev/null 2>&1; then
    echo "Error: neither pigz nor gzip found. One is required to compress output." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "Downloading ${ACCESSION} to $(pwd) ..."

# fasterq-dump produces uncompressed FASTQ files by default.
# --split-files writes R1 and R2 as separate files.
# --threads 4 is a reasonable default for most systems.
fasterq-dump \
    --split-files \
    --threads 4 \
    --progress \
    "${ACCESSION}"

echo "Compressing FASTQ files ..."

if command -v pigz >/dev/null 2>&1; then
    pigz -p 4 "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
else
    gzip "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
fi

echo ""
echo "Download complete."
echo ""
echo "Files:"
ls -lh "${ACCESSION}_1.fastq.gz" "${ACCESSION}_2.fastq.gz"
echo ""
echo "Notes:"
echo "  - UMI is in the read name, not in the first bases of the reads."
echo "  - Extract the 12 bp UMI with umi_tools extract before running kam."
echo "  - There is no SNV/indel truth set for this accession."
echo "  - See docs/benchmarking/05-public/umi_benchmark/README.md for preprocessing steps."
