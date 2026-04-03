#!/bin/sh
# Download SEQC2 HCC1395 somatic benchmark data.
#
# Downloads one Illumina WGS tumour-normal pair and the SNV/indel truth VCFs.
#
# Source: NCBI SRA (open access) + NCBI FTP (truth VCF).
# Paper: Xiao et al. "Achieving robust somatic mutation detection with deep
#        learning models derived from reference data sets of a cancer sample."
#        Genome Biology 22, 235 (2021).
#        https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02592-9
#
# Accessions downloaded:
#   SRR7890824  — HCC1395 tumour WGS  (~50 GB)
#   SRR7890827  — HCC1395BL normal WGS (~50 GB)
#   Truth VCFs  — ~50 MB
#
# Expected total disk space: ~100 GB.
# Runtime: 2–6 hours depending on network speed.
#
# Usage:
#   sh download_seqc2.sh [OUTPUT_DIR]
#
# If OUTPUT_DIR is not given, files are placed in ./data/seqc2/.

set -eu

TUMOUR_ACC="SRR7890824"
NORMAL_ACC="SRR7890827"
OUTPUT_DIR="${1:-data/seqc2}"

# NCBI FTP base for SEQC2 truth set files.
NCBI_FTP="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG"

# Check for required tools.
if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "Error: fasterq-dump not found." >&2
    echo "Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit" >&2
    exit 1
fi

if ! command -v wget >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1; then
    echo "Error: neither wget nor curl found. One is required to download truth VCFs." >&2
    exit 1
fi

if ! command -v pigz >/dev/null 2>&1 && ! command -v gzip >/dev/null 2>&1; then
    echo "Error: neither pigz nor gzip found. One is required to compress output." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Helper: download a URL with wget or curl.
download_url() {
    url="$1"
    outfile="$2"
    if command -v wget >/dev/null 2>&1; then
        wget -q --show-progress -O "${outfile}" "${url}"
    else
        curl -L --progress-bar -o "${outfile}" "${url}"
    fi
}

# Compress helper.
compress_fastq() {
    if command -v pigz >/dev/null 2>&1; then
        pigz -p 4 "$@"
    else
        gzip "$@"
    fi
}

echo "=== SEQC2 HCC1395 download ==="
echo "Output directory: $(pwd)"
echo ""

# --- Tumour ---
echo "Downloading tumour: ${TUMOUR_ACC} ..."
fasterq-dump \
    --split-files \
    --threads 4 \
    --progress \
    "${TUMOUR_ACC}"

echo "Compressing tumour FASTQ ..."
compress_fastq "${TUMOUR_ACC}_1.fastq" "${TUMOUR_ACC}_2.fastq"

# --- Normal ---
echo "Downloading normal: ${NORMAL_ACC} ..."
fasterq-dump \
    --split-files \
    --threads 4 \
    --progress \
    "${NORMAL_ACC}"

echo "Compressing normal FASTQ ..."
compress_fastq "${NORMAL_ACC}_1.fastq" "${NORMAL_ACC}_2.fastq"

# --- Truth VCFs ---
echo "Downloading truth VCFs from NCBI FTP ..."

# The short-form URLs (bit.ly) in the paper redirect to the NCBI FTP. Use the
# direct FTP path for stability; update if the NCBI FTP structure changes.
download_url \
    "${NCBI_FTP}/High-Confidence_Regions_v1.2/HighConf_SNVs_v1.2.vcf.gz" \
    "HCC1395_truth_SNV.vcf.gz"

download_url \
    "${NCBI_FTP}/High-Confidence_Regions_v1.2/HighConf_SNVs_v1.2.vcf.gz.tbi" \
    "HCC1395_truth_SNV.vcf.gz.tbi"

download_url \
    "${NCBI_FTP}/High-Confidence_Regions_v1.2/HighConf_Indels_v1.2.vcf.gz" \
    "HCC1395_truth_indel.vcf.gz"

download_url \
    "${NCBI_FTP}/High-Confidence_Regions_v1.2/HighConf_Indels_v1.2.vcf.gz.tbi" \
    "HCC1395_truth_indel.vcf.gz.tbi"

echo ""
echo "Download complete."
echo ""
echo "Files:"
ls -lh .
echo ""
echo "Notes:"
echo "  - SEQC2 data has no UMI. kam requires --no-umi mode (not yet implemented)."
echo "  - Evaluate against PASS variants in the truth VCF only."
echo "  - See docs/benchmarking/05-public/seqc2/README.md for evaluation steps."
echo ""
echo "  Truth VCF note: If the NCBI FTP paths above fail, check the current"
echo "  directory listing at:"
echo "  ${NCBI_FTP}/"
echo "  The original bit.ly short URLs from the paper are:"
echo "    SNV:   http://bit.ly/2DWuXzP"
echo "    Indel: http://bit.ly/2PVzzLn"
