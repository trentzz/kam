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

# Tracked config and metadata live in docs/.
ML_DIR="${REPO}/docs/project/experiments/01-ml-twist-duplex"

# Large simulation outputs live in bigdata/ (gitignored, mirrors Nextcloud).
BIGDATA_DIR="${REPO}/bigdata/experiments/01-ml-twist-duplex"
SIM_DIR="${BIGDATA_DIR}/simulations"
MODELS_DIR="${BIGDATA_DIR}/models"

# Load edit token from .env if present.
if [[ -f "${REPO}/.env" ]]; then
    # shellcheck source=/dev/null
    set -a
    source "${REPO}/.env"
    set +a
fi

if [[ -z "${NEXTCLOUD_EDIT_TOKEN:-}" ]]; then
    echo "[ERROR] NEXTCLOUD_EDIT_TOKEN is not set." >&2
    echo "        Copy .env.example to .env and fill in the token." >&2
    exit 1
fi

NEXTCLOUD_TOKEN="${NEXTCLOUD_EDIT_TOKEN}"
NEXTCLOUD_URL="https://nextcloudlocal.trentz.me/public.php/dav/files/${NEXTCLOUD_EDIT_TOKEN}"

mkdir -p ~/tmp

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
            local tarball="${HOME}/tmp/ml_twist_${split}_${batch_name}.tar.gz"

            echo "[TAR] ${split}/${batch_name} (${#batch_dirs[@]} dirs)"
            tar -czf "${tarball}" -C "${split_dir}" \
                --exclude="*_R1.fastq.gz" --exclude="*_R2.fastq.gz" \
                "${batch_dirs[@]}"
            upload_file "${tarball}" "experiments/01-ml-twist-duplex/${split}/${batch_name}.tar.gz"
            rm -f "${tarball}"

            batch_dirs=()
            (( batch_idx++ ))
        fi
    done < <(find "${split_dir}" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)

    # Upload any remaining directories
    if (( ${#batch_dirs[@]} > 0 )); then
        local batch_name
        batch_name="$(printf 'batch_%03d' "$batch_idx")"
        local tarball="${HOME}/tmp/ml_twist_${split}_${batch_name}.tar.gz"

        echo "[TAR] ${split}/${batch_name} (${#batch_dirs[@]} dirs)"
        tar -czf "${tarball}" -C "${split_dir}" \
            --exclude="*_R1.fastq.gz" --exclude="*_R2.fastq.gz" \
            "${batch_dirs[@]}"
        upload_file "${tarball}" "experiments/01-ml-twist-duplex/${split}/${batch_name}.tar.gz"
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
            upload_file "$f" "experiments/01-ml-twist-duplex/models/$(basename "$f")"
        done
    done
}

upload_features() {
    for csv in "${BIGDATA_DIR}"/*_features.csv.gz; do
        [[ -f "$csv" ]] || continue
        upload_file "$csv" "experiments/01-ml-twist-duplex/$(basename "$csv")"
    done
}

upload_configs() {
    local tarball="${HOME}/tmp/ml_twist_configs.tar.gz"
    echo "[TAR] configs (train + test)"
    tar -czf "${tarball}" -C "${REPO}" \
        docs/project/experiments/01-ml-twist-duplex/configs/
    upload_file "${tarball}" "experiments/01-ml-twist-duplex/configs.tar.gz"
    rm -f "${tarball}"
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
    configs)
        upload_configs
        ;;
    all)
        upload_configs
        upload_split_batches "train"
        upload_split_batches "test"
        upload_models
        upload_features
        ;;
    *)
        echo "[ERROR] Unknown scope '${SCOPE}'. Use: train | test | models | features | configs | all" >&2
        exit 1
        ;;
esac

echo ""
echo "=== Upload complete ==="
