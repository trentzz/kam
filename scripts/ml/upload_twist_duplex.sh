#!/usr/bin/env bash
# Upload Twist duplex ML dataset to Nextcloud.
#
# Uploads simulation batch tarballs, trained models, and feature CSVs.
#
# Usage:
#   ./scripts/ml/upload_twist_duplex.sh [train|test|models|features|all]
#
# Default: all
#
# Requires: curl

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ML_DIR="${REPO}/docs/benchmarking/ml-twist-duplex"
SIM_DIR="${ML_DIR}/simulations"
MODELS_DIR="${ML_DIR}/models"

NEXTCLOUD_URL="https://nextcloudlocal.trentz.me/public.php/dav/files/pTizAiSAJQsPcDo"
NEXTCLOUD_TOKEN="pTizAiSAJQsPcDo"

SCOPE="${1:-all}"

# ─── Helpers ─────────────────────────────────────────────────────────────────

upload_file() {
    local local_path="$1"
    local remote_path="$2"
    local filename
    filename="$(basename "$local_path")"

    echo "[UPLOAD] ${remote_path}"
    curl --silent --show-error \
        -T "${local_path}" \
        "${NEXTCLOUD_URL}/${remote_path}" \
        -u "${NEXTCLOUD_TOKEN}:"
    echo "  OK: ${filename}"
}

upload_split_batches() {
    local split="$1"
    local split_dir="${SIM_DIR}/${split}"

    if [[ ! -d "$split_dir" ]]; then
        echo "[WARN] No simulation directory for split '${split}': ${split_dir}" >&2
        return
    fi

    # Group sample directories into batches of 500 and upload
    local batch_idx=0
    local batch_dirs=()
    local count=0

    while IFS= read -r -d '' sample_dir; do
        batch_dirs+=("$(basename "$sample_dir")")
        (( count++ ))

        if (( count % 500 == 0 )); then
            local batch_name
            batch_name="$(printf 'batch_%03d' "$batch_idx")"
            local tarball="/tmp/ml_twist_${split}_${batch_name}.tar.gz"

            echo "[TAR] ${split}/${batch_name} (${#batch_dirs[@]} dirs)"
            tar -czf "${tarball}" -C "${split_dir}" "${batch_dirs[@]}"
            upload_file "${tarball}" "benchmarking/ml-twist-duplex/${split}/${batch_name}.tar.gz"
            rm -f "${tarball}"

            batch_dirs=()
            (( batch_idx++ ))
        fi
    done < <(find "${split_dir}" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)

    # Upload any remaining directories
    if (( ${#batch_dirs[@]} > 0 )); then
        local batch_name
        batch_name="$(printf 'batch_%03d' "$batch_idx")"
        local tarball="/tmp/ml_twist_${split}_${batch_name}.tar.gz"

        echo "[TAR] ${split}/${batch_name} (${#batch_dirs[@]} dirs)"
        tar -czf "${tarball}" -C "${split_dir}" "${batch_dirs[@]}"
        upload_file "${tarball}" "benchmarking/ml-twist-duplex/${split}/${batch_name}.tar.gz"
        rm -f "${tarball}"
    fi
}

upload_models() {
    if [[ ! -d "${MODELS_DIR}" ]]; then
        echo "[WARN] Models directory not found: ${MODELS_DIR}" >&2
        return
    fi
    for ext in txt json pkl onnx; do
        for f in "${MODELS_DIR}"/*.${ext}; do
            [[ -f "$f" ]] || continue
            upload_file "$f" "benchmarking/ml-twist-duplex/models/$(basename "$f")"
        done
    done
}

upload_features() {
    for csv in "${ML_DIR}"/*_features.csv.gz; do
        [[ -f "$csv" ]] || continue
        upload_file "$csv" "benchmarking/ml-twist-duplex/$(basename "$csv")"
    done
}

# ─── Main ────────────────────────────────────────────────────────────────────

echo "Nextcloud upload — Twist duplex ML dataset"
echo "Scope: ${SCOPE}"
echo ""

case "${SCOPE}" in
    train)
        upload_split_batches "train"
        ;;
    test)
        upload_split_batches "test"
        ;;
    models)
        upload_models
        ;;
    features)
        upload_features
        ;;
    all)
        upload_split_batches "train"
        upload_split_batches "test"
        upload_models
        upload_features
        ;;
    *)
        echo "[ERROR] Unknown scope '${SCOPE}'. Use: train | test | models | features | all" >&2
        exit 1
        ;;
esac

echo ""
echo "=== Upload complete ==="
