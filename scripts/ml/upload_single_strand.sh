#!/usr/bin/env bash
# Upload single-strand ML experiment data to Nextcloud.
#
# Uploads models, configs, samples, results, and training CSV for
# experiment 02-ml-single-strand.
#
# Usage:
#   ./scripts/ml/upload_single_strand.sh [models|configs|samples|results|training|all]
#
# Default: all
#
# Large files and tarballs are uploaded using Nextcloud's TUS-style chunked
# protocol (MKCOL → PUT chunks → MOVE) to stay under Cloudflare's 100 MB limit.
#
# Requires: curl, pigz (or gzip)

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIGDATA="${REPO}/bigdata/experiments/02-ml-single-strand"
NC_PREFIX="experiments/02-ml-single-strand"

# Load edit token from .env if present.
if [[ -f "${REPO}/.env" ]]; then
    set -a
    # shellcheck source=/dev/null
    source "${REPO}/.env"
    set +a
fi

if [[ -z "${NEXTCLOUD_EDIT_TOKEN:-}" ]]; then
    echo "[ERROR] NEXTCLOUD_EDIT_TOKEN is not set." >&2
    echo "        Copy .env.example to .env and fill in the token." >&2
    exit 1
fi

TOKEN="${NEXTCLOUD_EDIT_TOKEN}"
NC_URL="https://nextcloudlocal.trentz.me/public.php/dav/files/${TOKEN}"
NC_UPLOAD_URL="https://nextcloudlocal.trentz.me/public.php/dav/uploads/${TOKEN}"
CHUNK_SIZE=$((90 * 1024 * 1024))  # 90 MB chunks

mkdir -p ~/tmp

# ─── Helpers ─────────────────────────────────────────────────────────────────

# Ensure a remote directory exists, creating it and all parents if needed.
nc_mkdir_p() {
    local path="$1"  # e.g. "experiments/02-ml-single-strand/samples/test"
    local parts
    IFS='/' read -ra parts <<< "$path"
    local current=""
    for part in "${parts[@]}"; do
        [[ -z "$part" ]] && continue
        current="${current}${part}"
        local status
        status="$(curl -s -o /dev/null -w "%{http_code}" -X MKCOL \
            -u "${TOKEN}:" "${NC_URL}/${current}/")"
        # 201 = created, 405 = already exists — both are fine.
        if [[ "$status" != "201" && "$status" != "405" ]]; then
            echo "  [ERROR] MKCOL ${current}/ failed: HTTP ${status}" >&2
            return 1
        fi
        current="${current}/"
    done
}

# Upload a single file using chunked protocol.
chunked_upload() {
    local local_path="$1"
    local nc_dest="$2"  # path relative to share root
    local filename
    filename="$(basename "$local_path")"
    local file_size
    file_size="$(stat -c%s "$local_path")"
    local uuid
    uuid="ss-$(date +%s)-$$"
    local upload_base="${NC_UPLOAD_URL}/${uuid}"
    local dest_url="${NC_URL}/${nc_dest}"

    echo "[UPLOAD] ${nc_dest} ($(numfmt --to=iec "$file_size"))"

    # 1. Create upload slot.
    local status
    status="$(curl -s -o /dev/null -w "%{http_code}" -X MKCOL \
        -u "${TOKEN}:" "${upload_base}/")"
    if [[ "$status" != "201" ]]; then
        echo "  [ERROR] MKCOL failed: HTTP ${status}" >&2; return 1
    fi

    # 2. Upload chunks.
    local offset=0
    local chunk_idx=0
    while [[ $offset -lt $file_size ]]; do
        local chunk_path="${HOME}/tmp/ss_chunk_${uuid}_${chunk_idx}"
        dd if="$local_path" bs="${CHUNK_SIZE}" skip="$chunk_idx" count=1 \
            of="$chunk_path" 2>/dev/null
        local chunk_size
        chunk_size="$(stat -c%s "$chunk_path")"

        status="$(curl -s -o /dev/null -w "%{http_code}" \
            -T "$chunk_path" -u "${TOKEN}:" "${upload_base}/${offset}")"
        rm -f "$chunk_path"
        if [[ "$status" != "201" ]]; then
            echo "  [ERROR] chunk ${chunk_idx} failed: HTTP ${status}" >&2; return 1
        fi

        offset=$(( offset + chunk_size ))
        chunk_idx=$(( chunk_idx + 1 ))
        echo "  chunk ${chunk_idx} — $(numfmt --to=iec "$offset") / $(numfmt --to=iec "$file_size")"
    done

    # 3. Assemble.
    status="$(curl -s -o /dev/null -w "%{http_code}" -X MOVE \
        -u "${TOKEN}:" "${upload_base}/.file" \
        -H "Destination: ${dest_url}" \
        -H "Overwrite: T")"
    if [[ "$status" != "201" && "$status" != "204" ]]; then
        echo "  [ERROR] MOVE failed: HTTP ${status}" >&2; return 1
    fi
    echo "  OK: ${filename}"
}

# Tar a list of directories under a base dir and upload in chunks.
upload_dir_batches() {
    local base_dir="$1"       # local directory containing subdirs
    local nc_prefix="$2"      # nextcloud path prefix
    local batch_tag="$3"      # label for tarball names
    local batch_size="${4:-200}"

    if [[ ! -d "$base_dir" ]]; then
        echo "[WARN] Not found: ${base_dir}" >&2; return
    fi

    nc_mkdir_p "${nc_prefix}"

    local batch_idx=0
    local -a batch_dirs=()
    local count=0

    while IFS= read -r -d '' d; do
        batch_dirs+=("$(basename "$d")")
        count=$(( count + 1 ))
        if (( count % batch_size == 0 )); then
            local name
            name="$(printf 'batch_%03d' "$batch_idx")"
            local tarball="${HOME}/tmp/${batch_tag}_${name}.tar.gz"
            echo "[TAR] ${name} (${#batch_dirs[@]} dirs)"
            tar -czf "$tarball" -C "$base_dir" "${batch_dirs[@]}"
            chunked_upload "$tarball" "${nc_prefix}/${name}.tar.gz"
            rm -f "$tarball"
            batch_dirs=()
            batch_idx=$(( batch_idx + 1 ))
        fi
    done < <(find "$base_dir" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)

    if (( ${#batch_dirs[@]} > 0 )); then
        local name
        name="$(printf 'batch_%03d' "$batch_idx")"
        local tarball="${HOME}/tmp/${batch_tag}_${name}.tar.gz"
        echo "[TAR] ${name} (${#batch_dirs[@]} dirs)"
        tar -czf "$tarball" -C "$base_dir" "${batch_dirs[@]}"
        chunked_upload "$tarball" "${nc_prefix}/${name}.tar.gz"
        rm -f "$tarball"
    fi
}

upload_models() {
    local dir="${BIGDATA}/models"
    [[ -d "$dir" ]] || { echo "[WARN] No models dir" >&2; return; }
    nc_mkdir_p "${NC_PREFIX}/models"
    for f in "${dir}"/*.txt "${dir}"/*.json "${dir}"/*.pkl "${dir}"/*.onnx; do
        [[ -f "$f" ]] || continue
        chunked_upload "$f" "${NC_PREFIX}/models/$(basename "$f")"
    done
}

upload_configs() {
    local dir="${BIGDATA}/configs"
    [[ -d "$dir" ]] || { echo "[WARN] No configs dir" >&2; return; }
    nc_mkdir_p "${NC_PREFIX}"
    # Upload as a single tarball — configs are small YAML files.
    local tarball="${HOME}/tmp/ss_configs.tar.gz"
    echo "[TAR] configs"
    tar -czf "$tarball" -C "${BIGDATA}" configs/
    chunked_upload "$tarball" "${NC_PREFIX}/configs.tar.gz"
    rm -f "$tarball"
}

upload_samples() {
    # ML3 test/train splits live in subdirs.
    upload_dir_batches "${BIGDATA}/samples/test"  "${NC_PREFIX}/samples/test"  "ss_samples_test"  200
    upload_dir_batches "${BIGDATA}/samples/train" "${NC_PREFIX}/samples/train" "ss_samples_train" 200
    # ML2 samples are top-level dirs (exclude test/ and train/).
    local base="${BIGDATA}/samples"
    local nc_p="${NC_PREFIX}/samples"
    nc_mkdir_p "${nc_p}"
    local batch_idx=0 count=0
    local -a batch_dirs=()
    while IFS= read -r -d '' d; do
        local name
        name="$(basename "$d")"
        [[ "$name" == "test" || "$name" == "train" ]] && continue
        batch_dirs+=("$name")
        count=$(( count + 1 ))
        if (( count % 200 == 0 )); then
            local bname
            bname="$(printf 'batch_%03d' "$batch_idx")"
            local tarball="${HOME}/tmp/ss_samples_${bname}.tar.gz"
            echo "[TAR] samples/${bname} (${#batch_dirs[@]} dirs)"
            tar -czf "$tarball" -C "$base" "${batch_dirs[@]}"
            chunked_upload "$tarball" "${nc_p}/${bname}.tar.gz"
            rm -f "$tarball"
            batch_dirs=(); batch_idx=$(( batch_idx + 1 ))
        fi
    done < <(find "$base" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)
    if (( ${#batch_dirs[@]} > 0 )); then
        local bname
        bname="$(printf 'batch_%03d' "$batch_idx")"
        local tarball="${HOME}/tmp/ss_samples_${bname}.tar.gz"
        echo "[TAR] samples/${bname} (${#batch_dirs[@]} dirs)"
        tar -czf "$tarball" -C "$base" "${batch_dirs[@]}"
        chunked_upload "$tarball" "${nc_p}/${bname}.tar.gz"
        rm -f "$tarball"
    fi
}

upload_results() {
    # ML3 train/test subdirs.
    upload_dir_batches "${BIGDATA}/results/train" "${NC_PREFIX}/results/train" "ss_results_train" 200
    upload_dir_batches "${BIGDATA}/results/test"  "${NC_PREFIX}/results/test"  "ss_results_test"  200
    # Top-level ml2 result dirs (cam_* and sim_* at root of results/).
    local base="${BIGDATA}/results"
    local nc_p="${NC_PREFIX}/results"
    nc_mkdir_p "${nc_p}"
    local batch_idx=0 count=0
    local -a batch_dirs=()
    while IFS= read -r -d '' d; do
        local name
        name="$(basename "$d")"
        [[ "$name" == "test" || "$name" == "train" ]] && continue
        batch_dirs+=("$name")
        count=$(( count + 1 ))
        if (( count % 200 == 0 )); then
            local bname
            bname="$(printf 'batch_%03d' "$batch_idx")"
            local tarball="${HOME}/tmp/ss_results_${bname}.tar.gz"
            echo "[TAR] results/${bname} (${#batch_dirs[@]} dirs)"
            tar -czf "$tarball" -C "$base" "${batch_dirs[@]}"
            chunked_upload "$tarball" "${nc_p}/${bname}.tar.gz"
            rm -f "$tarball"
            batch_dirs=(); batch_idx=$(( batch_idx + 1 ))
        fi
    done < <(find "$base" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)
    if (( ${#batch_dirs[@]} > 0 )); then
        local bname
        bname="$(printf 'batch_%03d' "$batch_idx")"
        local tarball="${HOME}/tmp/ss_results_${bname}.tar.gz"
        echo "[TAR] results/${bname} (${#batch_dirs[@]} dirs)"
        tar -czf "$tarball" -C "$base" "${batch_dirs[@]}"
        chunked_upload "$tarball" "${nc_p}/${bname}.tar.gz"
        rm -f "$tarball"
    fi
}

upload_training() {
    local f="${BIGDATA}/training_data_v2.csv"
    [[ -f "$f" ]] || { echo "[WARN] training_data_v2.csv not found" >&2; return; }
    nc_mkdir_p "${NC_PREFIX}"
    chunked_upload "$f" "${NC_PREFIX}/training_data_v2.csv"
}

# ─── Main ────────────────────────────────────────────────────────────────────

SCOPE="${1:-all}"
echo "Nextcloud upload — single-strand ML (experiment 02)"
echo "Scope: ${SCOPE}"
echo ""

case "${SCOPE}" in
    models)   upload_models ;;
    configs)  upload_configs ;;
    samples)  upload_samples ;;
    results)  upload_results ;;
    training) upload_training ;;
    all)
        upload_models
        upload_configs
        upload_samples
        upload_results
        upload_training
        ;;
    *)
        echo "[ERROR] Unknown scope '${SCOPE}'. Use: models | configs | samples | results | training | all" >&2
        exit 1
        ;;
esac

echo ""
echo "=== Upload complete ==="
