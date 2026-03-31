"""Build ML training dataset v3 with an explicit train/test split.

Extends build_training_data_v2.py to cover the ML3 dataset structure:
  - docs/benchmarking/ml/samples/train/<sample>/   — ML3 training samples
  - docs/benchmarking/ml/samples/test/<sample>/    — ML3 test samples

The 'split' column is derived from the directory path: 'train' or 'test'
for ML3 samples, and 'legacy' for all earlier samples (snvindel + sv + ml).

Legacy samples are kept so the training set is as large as possible.
The test split contains only ML3 test samples to prevent contamination.

Outputs to docs/benchmarking/ml/training_data_v3.csv.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

# Legacy sample roots (all labelled split='legacy').
LEGACY_ROOTS = {
    "snvindel": REPO_ROOT / "docs/benchmarking/snvindel/samples",
    "sv":       REPO_ROOT / "docs/benchmarking/sv/samples",
    "ml":       REPO_ROOT / "docs/benchmarking/ml/samples",
}

# ML3 train and test roots.
ML3_TRAIN_ROOT = REPO_ROOT / "docs/benchmarking/ml/samples/train"
ML3_TEST_ROOT  = REPO_ROOT / "docs/benchmarking/ml/samples/test"

OUT_PATH = REPO_ROOT / "docs/benchmarking/ml/training_data_v3.csv"

TRUTH_KEY = ["chrom", "pos", "ref", "alt"]
CALL_COLS = [
    "chrom", "pos", "ref", "alt", "filter", "vaf", "vaf_lo", "vaf_hi",
    "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
]
PARAMS_COLS = [
    "params_coverage",
    "params_family_size_mean",
    "params_pcr_cycles",
    "params_fragment_mean",
]


def load_calls(path: Path) -> pd.DataFrame:
    """Load a calls TSV, returning an empty DataFrame on missing file."""
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    for col in CALL_COLS:
        if col not in df.columns:
            df[col] = float("nan")
    return df[CALL_COLS].copy()


def load_truth(path: Path) -> pd.DataFrame:
    """Load truth TSV. Requires chrom+pos+ref+alt columns."""
    if not path.exists():
        return pd.DataFrame(columns=TRUTH_KEY)
    df = pd.read_csv(path, sep="\t")
    for col in TRUTH_KEY:
        if col not in df.columns:
            raise ValueError(f"truth.tsv missing required column '{col}': {path}")
    return df[TRUTH_KEY].drop_duplicates()


def load_params(sample_dir: Path) -> dict:
    """Load params.json if present; return flattened numeric fields."""
    params_path = sample_dir / "params.json"
    if not params_path.exists():
        return {}
    with open(params_path) as f:
        raw = json.load(f)
    return {
        "params_coverage": raw.get("coverage", float("nan")),
        "params_family_size_mean": raw.get("family_size_mean", float("nan")),
        "params_pcr_cycles": raw.get("pcr_cycles", float("nan")),
        "params_fragment_mean": raw.get("fragment_mean", float("nan")),
        "split": raw.get("split", "legacy"),
    }


def derive_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add all derived feature columns and return the DataFrame."""
    nalt    = df["nalt"].fillna(0).astype(float)
    nref    = df["nref"].fillna(0).astype(float)
    ndupalt = df["ndupalt"].fillna(0).astype(float)
    nsimalt = df["nsimalt"].fillna(0).astype(float)
    vaf     = df["vaf"].fillna(0).astype(float)
    conf    = df["conf"].fillna(0).astype(float)
    sbp     = df["sbp"].fillna(1.0).astype(float)

    eps = 1e-9

    df["duplex_frac"]    = ndupalt / (nalt + eps)
    df["has_duplex"]     = (ndupalt > 0).astype(int)
    df["ci_width"]       = df["vaf_hi"] - df["vaf_lo"]
    df["ref_len"]        = df["ref"].fillna("").astype(str).str.len()
    df["alt_len"]        = df["alt"].fillna("").astype(str).str.len()
    df["alt_depth"]      = nalt + nref

    ref_len  = df["ref_len"].astype(float)
    alt_len  = df["alt_len"].astype(float)
    ci_width = df["ci_width"].fillna(0).astype(float)

    df["log_nalt"]       = np.log1p(nalt)
    df["log_nref"]       = np.log1p(nref)
    df["log_alt_depth"]  = np.log1p(nalt + nref)
    df["log_vaf"]        = np.log(vaf + 1e-6)

    df["vaf_times_conf"] = vaf * conf
    df["vaf_times_nalt"] = vaf * nalt
    df["nalt_over_conf"] = nalt / (conf + eps)
    df["ci_width_rel"]   = ci_width / (vaf + eps)
    df["snr"]            = nalt / (nref + 1.0)

    df["conf_sq"]        = conf ** 2
    df["nalt_sq"]        = nalt ** 2
    df["vaf_sq"]         = vaf ** 2

    df["ref_alt_len_ratio"] = ref_len / (alt_len + 1.0)
    df["indel_size"]     = (ref_len - alt_len).abs()

    df["duplex_enrichment"]  = ndupalt / (vaf * (nalt + nref) + eps)
    df["simplex_only_frac"]  = nsimalt / (nalt + eps)

    df["conf_above_99"]  = (conf > 0.99).astype(int)
    df["conf_above_999"] = (conf > 0.999).astype(int)
    df["sbp_above_05"]   = (sbp > 0.05).astype(int)

    # Encode variant class for the 33-feature rust-safe set.
    def classify_variant(row) -> str:
        r, a = row["ref_len"], row["alt_len"]
        if r == 1 and a == 1:
            return "SNV"
        mx = max(r, a)
        if mx < 50:
            return "Insertion" if a > r else "Deletion"
        return "LargeDeletion" if r > a else "TandemDuplication"

    vt_map = {
        "SNV": 0, "Insertion": 1, "Deletion": 2, "MNV": 3,
        "Complex": 4, "LargeDeletion": 5, "TandemDuplication": 6,
        "Inversion": 7, "Fusion": 8, "InvDel": 9, "NovelInsertion": 10,
    }
    df["variant_class"] = df.apply(classify_variant, axis=1)
    df["variant_class_enc"] = df["variant_class"].map(vt_map).fillna(0).astype(int)

    return df


def process_mode(
    sample_dir: Path,
    sample_id: str,
    dataset_type: str,
    mode: str,
    truth_df: pd.DataFrame,
    params: dict,
    split: str,
) -> pd.DataFrame:
    """Load one mode's TSV, label against truth, and return a tagged DataFrame."""
    tsv_name = "discovery.tsv" if mode == "discovery" else "tumour_informed.tsv"
    calls = load_calls(sample_dir / tsv_name)
    if calls.empty:
        return pd.DataFrame()

    truth_set = set(zip(
        truth_df["chrom"].astype(str),
        truth_df["pos"].astype(str),
        truth_df["ref"].astype(str),
        truth_df["alt"].astype(str),
    ))
    calls["label"] = calls.apply(
        lambda r: 1
        if (str(r["chrom"]), str(r["pos"]), str(r["ref"]), str(r["alt"])) in truth_set
        else 0,
        axis=1,
    )
    calls = derive_features(calls)
    calls["sample_id"]   = sample_id
    calls["dataset_type"] = dataset_type
    calls["mode"]        = mode
    calls["split"]       = split
    for col in PARAMS_COLS:
        calls[col] = params.get(col, float("nan"))
    return calls


def collect_samples(
    sample_root: Path,
    dataset_type: str,
    split: str,
    all_rows: list,
) -> None:
    """Process all sample directories under sample_root and append to all_rows."""
    if not sample_root.exists():
        return
    for sample_dir in sorted(p for p in sample_root.iterdir() if p.is_dir()):
        sample_id = sample_dir.name
        truth_df = load_truth(sample_dir / "truth.tsv")
        if truth_df.empty:
            continue
        params = load_params(sample_dir)
        # Use the split from params.json if present, else use the one passed in.
        effective_split = params.get("split", split)
        for mode in ("discovery", "tumour_informed"):
            df = process_mode(
                sample_dir, sample_id, dataset_type, mode,
                truth_df, params, effective_split,
            )
            if not df.empty:
                all_rows.append(df)


def main() -> None:
    all_rows: list[pd.DataFrame] = []

    # Legacy samples.
    for dataset_type, root in LEGACY_ROOTS.items():
        collect_samples(root, dataset_type, "legacy", all_rows)

    # ML3 train and test splits.
    collect_samples(ML3_TRAIN_ROOT, "ml3", "train", all_rows)
    collect_samples(ML3_TEST_ROOT,  "ml3", "test",  all_rows)

    if not all_rows:
        print("ERROR: no data collected", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_rows, ignore_index=True)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUT_PATH, index=False)
    print(f"Written {len(combined):,} rows to {OUT_PATH}")
    print(f"Columns: {len(combined.columns)}")

    print("\nSplit distribution:")
    split_stats = (
        combined.groupby(["split", "label"])
        .size()
        .rename("count")
        .reset_index()
    )
    print(split_stats.to_string(index=False))


if __name__ == "__main__":
    main()
