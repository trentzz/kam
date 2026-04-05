#!/usr/bin/env bash
# Time a conventional alignment pipeline on the same titration samples used by kam.
#
# Pipeline: BWA-MEM2 → samtools sort/index → GATK HaplotypeCaller (or
# freebayes as a lighter alternative).  Each stage is timed separately.
#
# This script is a TEMPLATE. The exact tool versions, reference genome path,
# and GATK call used in the comparison paper must be filled in before running.
# Sections that require user configuration are marked with TODO comments.
#
# Output:
#   docs/benchmarking/06-runtime/alignment_timings.csv
#
# CSV columns:
#   sample, wall_total_s, t_bwa_s, t_sort_s, t_index_s, t_gatk_s, peak_rss_mb
#
# Usage:
#   bash time_alignment.sh [--fastq-dir DIR] [--ref FILE] [--force]
#
#   --fastq-dir DIR   Directory containing per-sample FASTQ sub-directories.
#                     Defaults to $KAM_FASTQ_DIR or /data/titration-nondedup/fastqs.
#   --ref FILE        Reference genome FASTA (must be BWA-MEM2 indexed).
#   --force           Re-run samples that already have a timing entry.
#
# Prerequisites:
#   - bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2)
#   - samtools >= 1.17
#   - gatk >= 4.4 (or freebayes >= 1.3 as an alternative; see TODO below)
#   - /usr/bin/time (GNU time) for peak RSS measurement
#
# Notes:
#   - The reference genome is NOT included in this repository. Provide the
#     same reference that was used to design the capture panel.
#   - GATK HaplotypeCaller requires a .dict and .fai file alongside the FASTA.
#     Generate them with:
#       samtools faidx reference.fa
#       gatk CreateSequenceDictionary -R reference.fa
#   - Wall-clock time per stage is measured with 'date +%s%N'.
#   - GNU time (-v) records peak RSS from the kernel; this is more reliable
#     than self-reported metrics for multi-threaded tools.
#   - The pipeline uses 8 threads by default. Adjust THREADS as needed to
#     match the resources used for the kam timing run.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../" && pwd)"
EXPDIR="${REPO}/docs/benchmarking/runtime"
OUT_CSV="${EXPDIR}/alignment_timings.csv"

FASTQ_DIR="${KAM_FASTQ_DIR:-/data/titration-nondedup/fastqs}"
THREADS=8
FORCE=false

# TODO: Set this to the path of the indexed reference FASTA.
REF=""

# ── Argument parsing ──────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq-dir)  FASTQ_DIR="$2"; shift 2 ;;
        --ref)        REF="$2";       shift 2 ;;
        --force)      FORCE=true;     shift   ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ── Preflight checks ──────────────────────────────────────────────────────────

if [[ -z "$REF" ]]; then
    echo "[ERROR] --ref is required. Set it to the path of the BWA-MEM2-indexed reference FASTA." >&2
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference not found: $REF" >&2
    exit 1
fi

for tool in bwa-mem2 samtools gatk; do
    # TODO: Replace 'gatk' with 'freebayes' if using freebayes.
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[ERROR] Required tool not in PATH: $tool" >&2
        exit 1
    fi
done

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "[ERROR] FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

GNU_TIME=""
if /usr/bin/time --version >/dev/null 2>&1; then
    GNU_TIME="/usr/bin/time"
fi

# ── CSV header ────────────────────────────────────────────────────────────────

if [[ ! -f "$OUT_CSV" ]]; then
    echo "sample,wall_total_s,t_bwa_s,t_sort_s,t_index_s,t_gatk_s,peak_rss_mb" \
        > "$OUT_CSV"
fi

# ── Timing helper ──────────────────────────────────────────────────────────────

# Usage: elapsed_s START_NS
# Returns elapsed seconds (3 decimal places) since START_NS.
elapsed_s() {
    local start="$1"
    local now
    now=$(date +%s%N)
    echo "scale=3; ($now - $start) / 1000000000" | bc
}

# ── Main loop ─────────────────────────────────────────────────────────────────

FAILED=()

for sample_dir in "$FASTQ_DIR"/*/; do
    [[ -d "$sample_dir" ]] || continue
    sample="$(basename "$sample_dir")"

    if [[ "$FORCE" == "false" ]] && grep -q "^${sample}," "$OUT_CSV" 2>/dev/null; then
        echo "[SKIP] ${sample} (already in CSV)"
        continue
    fi

    r1=$(ls "${sample_dir}"*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sample_dir}"*_R2.fastq.gz 2>/dev/null | head -1)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[WARN] No FASTQ pair in ${sample_dir}" >&2
        FAILED+=("$sample")
        continue
    fi

    work_dir="${EXPDIR}/tmp_align_${sample}"
    rm -rf "$work_dir"
    mkdir -p "$work_dir"

    echo "[ALIGN] ${sample}"

    t_total_start=$(date +%s%N)
    peak_rss_mb=""

    # ── Stage 1: BWA-MEM2 alignment ──────────────────────────────────────────
    t_bwa_start=$(date +%s%N)

    bam_unsorted="${work_dir}/unsorted.bam"
    TIME_OUT="${work_dir}/time_bwa.txt"

    if [[ -n "$GNU_TIME" ]]; then
        if ! "$GNU_TIME" -v -o "$TIME_OUT" \
                bwa-mem2 mem -t "$THREADS" \
                    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
                    "$REF" "$r1" "$r2" \
                | samtools view -bS -o "$bam_unsorted" 2>/dev/null; then
            echo "[ERROR] BWA-MEM2 failed for ${sample}" >&2
            FAILED+=("$sample")
            rm -rf "$work_dir"
            continue
        fi
        peak_rss_mb=$(grep "Maximum resident set size" "$TIME_OUT" \
            | awk '{print $NF / 1024}' 2>/dev/null || echo "")
    else
        if ! bwa-mem2 mem -t "$THREADS" \
                -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
                "$REF" "$r1" "$r2" \
            | samtools view -bS -o "$bam_unsorted" 2>/dev/null; then
            echo "[ERROR] BWA-MEM2 failed for ${sample}" >&2
            FAILED+=("$sample")
            rm -rf "$work_dir"
            continue
        fi
    fi

    t_bwa_s=$(elapsed_s "$t_bwa_start")

    # ── Stage 2: samtools sort ────────────────────────────────────────────────
    t_sort_start=$(date +%s%N)
    bam_sorted="${work_dir}/sorted.bam"

    if ! samtools sort -@ "$THREADS" -o "$bam_sorted" "$bam_unsorted" 2>/dev/null; then
        echo "[ERROR] samtools sort failed for ${sample}" >&2
        FAILED+=("$sample")
        rm -rf "$work_dir"
        continue
    fi
    rm -f "$bam_unsorted"
    t_sort_s=$(elapsed_s "$t_sort_start")

    # ── Stage 3: samtools index ───────────────────────────────────────────────
    t_index_start=$(date +%s%N)

    if ! samtools index "$bam_sorted" 2>/dev/null; then
        echo "[ERROR] samtools index failed for ${sample}" >&2
        FAILED+=("$sample")
        rm -rf "$work_dir"
        continue
    fi
    t_index_s=$(elapsed_s "$t_index_start")

    # ── Stage 4: GATK HaplotypeCaller ────────────────────────────────────────
    # TODO: Adjust the GATK command to match the actual calling mode used.
    # For somatic calling (Mutect2), replace HaplotypeCaller with Mutect2 and
    # add the appropriate germline resource and panel-of-normals arguments.
    t_gatk_start=$(date +%s%N)
    vcf_out="${work_dir}/calls.vcf.gz"

    if ! gatk HaplotypeCaller \
            -R "$REF" \
            -I "$bam_sorted" \
            -O "$vcf_out" \
            --native-pair-hmm-threads "$THREADS" \
            2>/dev/null; then
        echo "[ERROR] GATK failed for ${sample}" >&2
        FAILED+=("$sample")
        rm -rf "$work_dir"
        continue
    fi
    t_gatk_s=$(elapsed_s "$t_gatk_start")

    # ── Total wall time ───────────────────────────────────────────────────────
    wall_total_s=$(elapsed_s "$t_total_start")

    echo "${sample},${wall_total_s},${t_bwa_s},${t_sort_s},${t_index_s},${t_gatk_s},${peak_rss_mb}" \
        >> "$OUT_CSV"

    rm -rf "$work_dir"
done

echo ""
echo "=== DONE ==="
echo "Output: ${OUT_CSV}"
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
