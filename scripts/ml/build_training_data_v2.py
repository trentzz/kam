"""Build enhanced ML training dataset with richer derived features.

Extends build_training_data.py with log transforms, ratios, interactions,
binned features, and confidence flags. Also ingests params.json from each
sample directory when present.

Outputs to bigdata/experiments/02-ml-single-strand/training_data_v2.csv.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Repo root is two levels up from this script (scripts/ml/)
REPO_ROOT = Path(__file__).resolve().parents[2]

SAMPLE_ROOTS = {
    "snvindel": REPO_ROOT / "docs/benchmarking/01-snvindel/samples",
    "sv": REPO_ROOT / "docs/benchmarking/02-sv-core/samples",
    "ml": REPO_ROOT / "docs/project/experiments/02-ml-single-strand/samples",
}

OUT_PATH = REPO_ROOT / "bigdata/experiments/02-ml-single-strand/training_data_v2.csv"

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
    """Load truth TSV. Require chrom+pos+ref+alt columns."""
    if not path.exists():
        return pd.DataFrame(columns=TRUTH_KEY)
    df = pd.read_csv(path, sep="\t")
    for col in TRUTH_KEY:
        if col not in df.columns:
            raise ValueError(f"truth.tsv missing required column '{col}': {path}")
    return df[TRUTH_KEY].drop_duplicates()


def load_params(sample_dir: Path) -> dict:
    """Load params.json if present; return dict of flattened param values."""
    params_path = sample_dir / "params.json"
    if not params_path.exists():
        return {}
    with open(params_path) as f:
        raw = json.load(f)
    # Flatten top-level numeric keys we care about
    keys = {
        "coverage": "params_coverage",
        "family_size_mean": "params_family_size_mean",
        "pcr_cycles": "params_pcr_cycles",
        "fragment_mean": "params_fragment_mean",
    }
    return {v: raw.get(k, float("nan")) for k, v in keys.items()}


def derive_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add all derived feature columns and return the DataFrame."""
    nalt = df["nalt"].fillna(0).astype(float)
    nref = df["nref"].fillna(0).astype(float)
    ndupalt = df["ndupalt"].fillna(0).astype(float)
    nsimalt = df["nsimalt"].fillna(0).astype(float)
    vaf = df["vaf"].fillna(0).astype(float)
    conf = df["conf"].fillna(0).astype(float)
    sbp = df["sbp"].fillna(1.0).astype(float)

    # ---- Base features ----
    df["duplex_frac"] = ndupalt.where(nalt > 0, 0.0) / nalt.where(nalt > 0, 1.0)
    df["has_duplex"] = (ndupalt > 0).astype(int)
    df["ci_width"] = df["vaf_hi"] - df["vaf_lo"]
    df["ref_len"] = df["ref"].fillna("").astype(str).str.len()
    df["alt_len"] = df["alt"].fillna("").astype(str).str.len()
    df["alt_depth"] = nalt + nref

    ref_len = df["ref_len"].astype(float)
    alt_len = df["alt_len"].astype(float)
    max_len = df[["ref_len", "alt_len"]].max(axis=1)
    ci_width = df["ci_width"].fillna(0).astype(float)

    def classify(row):
        if row["ref_len"] == 1 and row["alt_len"] == 1:
            return "SNV"
        if max_len[row.name] < 50:
            return "indel"
        return "SV"

    df["variant_class"] = df.apply(classify, axis=1)

    # ---- Log transforms ----
    df["log_nalt"] = np.log1p(nalt)
    df["log_nref"] = np.log1p(nref)
    df["log_alt_depth"] = np.log1p(nalt + nref)
    df["log_vaf"] = np.log(vaf + 1e-6)

    # ---- Ratios and interactions ----
    raw_alt_frac = nalt / (nalt + nref + 1e-9)
    df["observed_vaf_diff"] = vaf - raw_alt_frac
    df["vaf_times_conf"] = vaf * conf
    df["vaf_times_nalt"] = vaf * nalt
    df["nalt_over_conf"] = nalt / (conf + 1e-9)
    df["ci_width_rel"] = ci_width / (vaf + 1e-9)
    df["snr"] = nalt / (nref + 1.0)

    # ---- Squared terms ----
    df["conf_sq"] = conf ** 2
    df["nalt_sq"] = nalt ** 2
    df["vaf_sq"] = vaf ** 2

    # ---- Length-based features ----
    df["ref_alt_len_ratio"] = ref_len / (alt_len + 1e-9)
    df["indel_size"] = (ref_len - alt_len).abs()

    # ---- Binned features ----
    alt_depth = df["alt_depth"].astype(float)
    # depth_bucket: 5 quantile bins (0–4)
    try:
        df["depth_bucket"] = pd.qcut(
            alt_depth, q=5, labels=False, duplicates="drop"
        ).fillna(0).astype(int)
    except ValueError:
        df["depth_bucket"] = 0

    # sbp_cat: <0.01 → 0, 0.01–0.5 → 1, >=0.5 → 2
    sbp_cat = np.where(sbp < 0.01, 0, np.where(sbp < 0.5, 1, 2))
    df["sbp_cat"] = sbp_cat.astype(int)

    # ---- Duplex enrichment features ----
    df["duplex_enrichment"] = ndupalt / (nalt + 1e-9)
    df["simplex_only_frac"] = nsimalt / (nalt + 1e-9)

    # ---- Confidence threshold flags ----
    df["conf_above_99"] = (conf > 0.99).astype(int)
    df["conf_above_999"] = (conf > 0.999).astype(int)
    df["sbp_above_05"] = (sbp > 0.05).astype(int)

    return df


def process_mode(
    sample_dir: Path,
    sample_id: str,
    dataset_type: str,
    mode: str,
    truth_df: pd.DataFrame,
    params: dict,
) -> pd.DataFrame:
    """Load one mode's TSV, label against truth, and return a tagged DataFrame."""
    tsv_name = "discovery.tsv" if mode == "discovery" else "tumour_informed.tsv"
    calls = load_calls(sample_dir / tsv_name)

    if calls.empty:
        return pd.DataFrame()

    truth_set = set(
        zip(
            truth_df["chrom"].astype(str),
            truth_df["pos"].astype(str),
            truth_df["ref"].astype(str),
            truth_df["alt"].astype(str),
        )
    )

    calls["label"] = calls.apply(
        lambda r: 1
        if (str(r["chrom"]), str(r["pos"]), str(r["ref"]), str(r["alt"])) in truth_set
        else 0,
        axis=1,
    )

    calls = derive_features(calls)
    calls["sample_id"] = sample_id
    calls["dataset_type"] = dataset_type
    calls["mode"] = mode

    # Params columns — fill NaN if not available
    for col in PARAMS_COLS:
        calls[col] = params.get(col, float("nan"))

    return calls


def main():
    all_rows = []

    for dataset_type, sample_root in SAMPLE_ROOTS.items():
        if not sample_root.exists():
            print(f"WARNING: {sample_root} does not exist — skipping", file=sys.stderr)
            continue

        sample_dirs = sorted(p for p in sample_root.iterdir() if p.is_dir())

        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            truth_df = load_truth(sample_dir / "truth.tsv")

            if truth_df.empty:
                print(f"WARNING: no truth for {sample_id} — skipping", file=sys.stderr)
                continue

            params = load_params(sample_dir)

            for mode in ("discovery", "tumour_informed"):
                df = process_mode(sample_dir, sample_id, dataset_type, mode, truth_df, params)
                if not df.empty:
                    all_rows.append(df)

    if not all_rows:
        print("ERROR: no data collected", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_rows, ignore_index=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUT_PATH, index=False)
    print(f"Written {len(combined):,} rows to {OUT_PATH}")
    print(f"Columns: {len(combined.columns)}")

    # Label statistics
    print("\nLabel statistics:")
    stats = (
        combined.groupby(["dataset_type", "mode", "label"])
        .size()
        .rename("count")
        .reset_index()
    )
    print(stats.to_string(index=False))

    total_pos = (combined["label"] == 1).sum()
    total_neg = (combined["label"] == 0).sum()
    total = len(combined)
    print(f"\nTotal rows:                      {total:,}")
    print(f"Total positives (true variants): {total_pos:,}")
    print(f"Total negatives (false calls):   {total_neg:,}")
    print(f"Positive rate:                   {total_pos / total:.5f} ({total_pos / total * 100:.4f}%)")

    print("\nNew feature columns (sample):")
    new_cols = [
        c for c in combined.columns
        if c not in [
            "chrom", "pos", "ref", "alt", "filter", "vaf", "vaf_lo", "vaf_hi",
            "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
            "label", "sample_id", "dataset_type", "mode",
            "duplex_frac", "has_duplex", "ci_width", "ref_len", "alt_len",
            "alt_depth", "variant_class",
        ]
    ]
    print(", ".join(new_cols))


if __name__ == "__main__":
    main()
