#!/bin/sh
# Download COLO829 Illumina WGS data and the somatic SV truth set.
#
# Downloads the Illumina HiSeq X Ten tumour-normal pair from ENA and the
# 68-variant SV truth VCF from Zenodo.
#
# Source:
#   Raw FASTQ:  ENA project PRJEB27698 (open access, no restrictions).
#   Truth VCF:  Zenodo record 4716169.
#
# Paper: van Dijk et al. "A multi-platform reference for somatic structural
#        variation detection." Cell Genomics 2, 100082 (2022).
#        https://www.sciencedirect.com/science/article/pii/S2666979X22000726
#
# Accessions downloaded:
#   ERR1341793  — COLO829 tumour Illumina WGS  (~100 GB)
#   ERR1341794  — COLO829BL normal Illumina WGS (~80 GB)
#   Truth VCFs  — <1 MB
#
# Expected total disk space: ~180 GB.
# Runtime: 3–8 hours depending on network speed.
#
# Usage:
#   sh download_colo829.sh [OUTPUT_DIR]
#
# If OUTPUT_DIR is not given, files are placed in ./data/colo829/.

set -eu

TUMOUR_ACC="ERR1341793"
NORMAL_ACC="ERR1341794"
OUTPUT_DIR="${1:-data/colo829}"

# Zenodo record for the COLO829 truth set.
ZENODO_BASE="https://zenodo.org/records/4716169/files"

# ENA FTP base for paired-end FASTQ files.
# ENA stores files at: ftp.sra.ebi.ac.uk/vol1/fastq/<6-char-prefix>/<accession>/
ena_ftp_url() {
    acc="$1"
    file="$2"
    # ENA organises by first 6 chars of accession number (drop the ERR prefix).
    prefix="$(echo "${acc}" | cut -c 1-6)"
    echo "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00$(echo "${acc}" | tail -c 2)/${acc}/${file}"
}

# Check for required tools.
if ! command -v wget >/dev/null 2>&1 && ! command -v curl >/dev/null 2>&1; then
    echo "Error: neither wget nor curl found." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Helper: download a URL with wget or curl.
download_url() {
    url="$1"
    outfile="$2"
    echo "  ${url} -> ${outfile}"
    if command -v wget >/dev/null 2>&1; then
        wget -q --show-progress -O "${outfile}" "${url}"
    else
        curl -L --progress-bar -o "${outfile}" "${url}"
    fi
}

echo "=== COLO829 somatic SV benchmark download ==="
echo "Output directory: $(pwd)"
echo ""

# --- Tumour Illumina WGS ---
echo "Downloading tumour: ${TUMOUR_ACC} ..."
download_url \
    "$(ena_ftp_url "${TUMOUR_ACC}" "${TUMOUR_ACC}_1.fastq.gz")" \
    "${TUMOUR_ACC}_1.fastq.gz"
download_url \
    "$(ena_ftp_url "${TUMOUR_ACC}" "${TUMOUR_ACC}_2.fastq.gz")" \
    "${TUMOUR_ACC}_2.fastq.gz"

# --- Normal Illumina WGS ---
echo "Downloading normal: ${NORMAL_ACC} ..."
download_url \
    "$(ena_ftp_url "${NORMAL_ACC}" "${NORMAL_ACC}_1.fastq.gz")" \
    "${NORMAL_ACC}_1.fastq.gz"
download_url \
    "$(ena_ftp_url "${NORMAL_ACC}" "${NORMAL_ACC}_2.fastq.gz")" \
    "${NORMAL_ACC}_2.fastq.gz"

# --- Truth VCFs from Zenodo ---
echo "Downloading truth VCFs from Zenodo record 4716169 ..."

# GRCh37 truth VCF (primary).
download_url \
    "${ZENODO_BASE}/COLO829_truthset_somatic_v4.1.vcf?download=1" \
    "COLO829_truthset_somatic_v4.1.vcf"

# hg38 liftover truth VCF.
download_url \
    "${ZENODO_BASE}/COLO829_truthset_somatic_hg38.vcf?download=1" \
    "COLO829_truthset_somatic_hg38.vcf"

echo ""
echo "Download complete."
echo ""
echo "Files:"
ls -lh .
echo ""
echo "Truth set summary: 68 validated somatic SVs"
echo "  38 DEL, 13 TRA, 7 DUP, 7 INV, 3 INS"
echo ""
echo "Notes:"
echo "  - COLO829 data has no UMI. kam requires --no-umi mode (not yet implemented)."
echo "  - Use GRCh37 truth VCF unless you align to hg38."
echo "  - Evaluate with Truvari: see docs/benchmarking/05-public/colo829/README.md."
echo "  - Confirm ENA accessions are current at:"
echo "    https://www.ebi.ac.uk/ena/browser/view/PRJEB27698"
echo ""
echo "  Zenodo record: https://zenodo.org/records/4716169"
echo "  GitHub:        https://github.com/UMCUGenetics/COLO829_somaticSV"
