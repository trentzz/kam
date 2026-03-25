#!/usr/bin/env bash
# Time kam on a representative set of titration samples.
#
# Runs kam in discovery mode on each sample, records wall-clock time, and
# extracts per-stage timings from the QC JSON emitted by kam.
#
# Output:
#   docs/benchmarking/runtime/kam_timings.csv
#
# CSV columns:
#   sample, wall_time_s, t_assemble_ms, t_index_ms, t_pathfind_ms,
#   t_call_ms, t_output_ms, peak_rss_mb
#
# Usage:
#   bash time_kam.sh [--fastq-dir DIR] [--targets FILE] [--force]
#
#   --fastq-dir DIR   Directory containing per-sample FASTQ sub-directories.
#                     Defaults to $KAM_FASTQ_DIR or /data/titration-nondedup/fastqs.
#   --targets FILE    Targets FASTA passed to kam. Defaults to
#                     docs/benchmarking/snvindel/scripts/targets_100bp.fa.
#   --force           Re-run samples that already have a timing entry.
#
# Prerequisites:
#   - target/release/kam must be built
#   - /usr/bin/time (GNU time) or the shell built-in 'time' with %e support
#   - python3 (for QC JSON parsing)
#
# Notes:
#   - The script uses the GNU 'time' command (/usr/bin/time -v) to measure
#     peak RSS. If unavailable, peak_rss_mb is left blank.
#   - Per-stage timings are read from the qc.json file written by kam to the
#     output directory. The expected JSON path is:
#       <out_dir>/qc.json
#     with fields: stage_timings.assemble_ms, stage_timings.index_ms,
#     stage_timings.pathfind_ms, stage_timings.call_ms, stage_timings.output_ms
#   - Wall-clock time is measured with 'date +%s%N' for sub-second precision.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../" && pwd)"
KAM="${REPO}/target/release/kam"
EXPDIR="${REPO}/docs/benchmarking/runtime"
OUT_CSV="${EXPDIR}/kam_timings.csv"

FASTQ_DIR="${KAM_FASTQ_DIR:-/data/titration-nondedup/fastqs}"
TARGETS="${REPO}/docs/benchmarking/snvindel/scripts/targets_100bp.fa"
FORCE=false

# ── Argument parsing ──────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq-dir)  FASTQ_DIR="$2"; shift 2 ;;
        --targets)    TARGETS="$2"; shift 2 ;;
        --force)      FORCE=true; shift ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# ── Preflight checks ──────────────────────────────────────────────────────────

if [[ ! -x "$KAM" ]]; then
    echo "[ERROR] kam binary not found: $KAM" >&2
    echo "        Build with: cargo build --release" >&2
    exit 1
fi

if [[ ! -f "$TARGETS" ]]; then
    echo "[ERROR] Targets file not found: $TARGETS" >&2
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "[ERROR] FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

# Check for GNU time (provides peak RSS).
GNU_TIME=""
if /usr/bin/time --version >/dev/null 2>&1; then
    GNU_TIME="/usr/bin/time"
fi

# ── CSV header ────────────────────────────────────────────────────────────────

if [[ ! -f "$OUT_CSV" ]]; then
    echo "sample,wall_time_s,t_assemble_ms,t_index_ms,t_pathfind_ms,t_call_ms,t_output_ms,peak_rss_mb" \
        > "$OUT_CSV"
fi

# ── Python helper: extract QC JSON fields ─────────────────────────────────────
# Written to a temp file to avoid a subprocess per sample.

QC_PARSER=$(mktemp /tmp/parse_qc_XXXXXX.py)
trap 'rm -f "$QC_PARSER"' EXIT

cat > "$QC_PARSER" << 'PYEOF'
#!/usr/bin/env python3
"""Read a kam qc.json and print stage timings as CSV fields.

Usage: python3 parse_qc.py <qc_json>

Output (one line, no header):
  t_assemble_ms,t_index_ms,t_pathfind_ms,t_call_ms,t_output_ms,peak_rss_mb
"""
import json
import sys

path = sys.argv[1]
try:
    with open(path) as fh:
        data = json.load(fh)
except (FileNotFoundError, json.JSONDecodeError) as e:
    print(f",,,,,,", file=sys.stdout)
    sys.exit(0)

st = data.get("stage_timings", {})
fields = [
    st.get("assemble_ms", ""),
    st.get("index_ms",    ""),
    st.get("pathfind_ms", ""),
    st.get("call_ms",     ""),
    st.get("output_ms",   ""),
    data.get("peak_rss_mb", ""),
]
print(",".join(str(f) for f in fields))
PYEOF

# ── Main loop ─────────────────────────────────────────────────────────────────

FAILED=()

# Find all sample directories (each should contain *_R1.fastq.gz and *_R2.fastq.gz).
for sample_dir in "$FASTQ_DIR"/*/; do
    [[ -d "$sample_dir" ]] || continue
    sample="$(basename "$sample_dir")"

    # Skip if already timed and not forcing a re-run.
    if [[ "$FORCE" == "false" ]] && grep -q "^${sample}," "$OUT_CSV" 2>/dev/null; then
        echo "[SKIP] ${sample} (already in CSV; use --force to re-run)"
        continue
    fi

    r1=$(ls "${sample_dir}"*_R1.fastq.gz 2>/dev/null | head -1)
    r2=$(ls "${sample_dir}"*_R2.fastq.gz 2>/dev/null | head -1)

    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[WARN] No FASTQ pair found in ${sample_dir}" >&2
        FAILED+=("$sample")
        continue
    fi

    out_dir="${EXPDIR}/tmp_timing_${sample}"
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    echo "[RUN] ${sample}"

    # Record start time (nanoseconds for precision).
    t_start=$(date +%s%N)

    if [[ -n "$GNU_TIME" ]]; then
        # GNU time writes to stderr; capture it.
        TIME_OUT=$(mktemp /tmp/gnu_time_XXXXXX.txt)
        if ! "$GNU_TIME" -v -o "$TIME_OUT" \
                "$KAM" run \
                    --r1 "$r1" --r2 "$r2" \
                    --targets "$TARGETS" \
                    --output-dir "$out_dir" \
                    --output-format vcf,tsv \
                    2>/dev/null; then
            echo "[ERROR] kam failed for ${sample}" >&2
            FAILED+=("$sample")
            rm -rf "$out_dir" "$TIME_OUT"
            continue
        fi
        peak_rss=$(grep "Maximum resident set size" "$TIME_OUT" \
            | awk '{print $NF / 1024}' 2>/dev/null || echo "")
        rm -f "$TIME_OUT"
    else
        if ! "$KAM" run \
                --r1 "$r1" --r2 "$r2" \
                --targets "$TARGETS" \
                --output-dir "$out_dir" \
                --output-format vcf,tsv \
                2>/dev/null; then
            echo "[ERROR] kam failed for ${sample}" >&2
            FAILED+=("$sample")
            rm -rf "$out_dir"
            continue
        fi
        peak_rss=""
    fi

    t_end=$(date +%s%N)
    wall_time_s=$(echo "scale=3; ($t_end - $t_start) / 1000000000" | bc)

    # Extract per-stage timings from the QC JSON.
    qc_json="${out_dir}/qc.json"
    stage_fields=$(python3 "$QC_PARSER" "$qc_json" 2>/dev/null || echo ",,,,,")

    # If GNU time already gave us RSS, prefer it (more accurate than QC JSON).
    if [[ -n "$peak_rss" ]]; then
        # Replace the last field with the GNU time value.
        stage_fields=$(echo "$stage_fields" | awk -F',' -v rss="$peak_rss" \
            '{OFS=","; $NF=rss; print}')
    fi

    echo "${sample},${wall_time_s},${stage_fields}" >> "$OUT_CSV"

    rm -rf "$out_dir"
done

echo ""
echo "=== DONE ==="
echo "Output: ${OUT_CSV}"
echo "Failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "  $f"; done
