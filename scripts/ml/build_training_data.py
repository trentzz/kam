"""Build ML training dataset from benchmarking sample directories.

Walks snvindel and sv sample directories, loads discovery.tsv,
tumour_informed.tsv, and truth.tsv for each sample, matches calls
to truth by exact chrom+pos+ref+alt, derives features, and writes
a combined CSV to bigdata/experiments/02-ml-single-strand/training_data.csv.
"""

import os
import sys
from pathlib import Path

import pandas as pd

# Repo root is two levels up from this script (scripts/ml/)
REPO_ROOT = Path(__file__).resolve().parents[2]

SAMPLE_ROOTS = {
    "snvindel": REPO_ROOT / "docs/benchmarking/01-snvindel/samples",
    "sv": REPO_ROOT / "docs/benchmarking/02-sv-core/samples",
}

OUT_PATH = REPO_ROOT / "bigdata/experiments/02-ml-single-strand/training_data.csv"

TRUTH_KEY = ["chrom", "pos", "ref", "alt"]
CALL_COLS = ["chrom", "pos", "ref", "alt", "filter", "vaf", "vaf_lo", "vaf_hi",
             "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf"]


def load_calls(path: Path) -> pd.DataFrame:
    """Load a calls TSV, returning an empty DataFrame on missing file."""
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    # Ensure expected columns exist; fill missing with NaN
    for col in CALL_COLS:
        if col not in df.columns:
            df[col] = float("nan")
    return df[CALL_COLS]


def load_truth(path: Path) -> pd.DataFrame:
    """Load truth TSV with at minimum chrom+pos+ref+alt columns."""
    if not path.exists():
        return pd.DataFrame(columns=TRUTH_KEY)
    df = pd.read_csv(path, sep="\t")
    for col in TRUTH_KEY:
        if col not in df.columns:
            raise ValueError(f"truth.tsv missing required column '{col}': {path}")
    return df[TRUTH_KEY].drop_duplicates()


def derive_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add derived feature columns in-place and return the DataFrame."""
    nalt = df["nalt"].fillna(0)
    ndupalt = df["ndupalt"].fillna(0)

    df["duplex_frac"] = ndupalt.where(nalt > 0, 0) / nalt.where(nalt > 0, 1)
    df["has_duplex"] = (ndupalt > 0).astype(int)
    df["ci_width"] = df["vaf_hi"] - df["vaf_lo"]
    df["ref_len"] = df["ref"].str.len()
    df["alt_len"] = df["alt"].str.len()
    df["alt_depth"] = df["nalt"].fillna(0) + df["nref"].fillna(0)

    max_len = df[["ref_len", "alt_len"]].max(axis=1)

    def classify(row):
        if row["ref_len"] == 1 and row["alt_len"] == 1:
            return "SNV"
        if max_len[row.name] < 50:
            return "indel"
        return "SV"

    df["variant_class"] = df.apply(classify, axis=1)
    return df


def process_mode(
    sample_dir: Path,
    sample_id: str,
    dataset_type: str,
    mode: str,
    truth_df: pd.DataFrame,
) -> pd.DataFrame:
    """Load one mode's TSV, label against truth, and return a tagged DataFrame."""
    tsv_name = "discovery.tsv" if mode == "discovery" else "tumour_informed.tsv"
    calls = load_calls(sample_dir / tsv_name)

    if calls.empty:
        return pd.DataFrame()

    # Match calls to truth
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

            for mode in ("discovery", "tumour_informed"):
                df = process_mode(sample_dir, sample_id, dataset_type, mode, truth_df)
                if not df.empty:
                    all_rows.append(df)

    if not all_rows:
        print("ERROR: no data collected", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_rows, ignore_index=True)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUT_PATH, index=False)
    print(f"Written {len(combined):,} rows to {OUT_PATH}")

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
    print(f"\nTotal positives (true variants): {total_pos:,}")
    print(f"Total negatives (false calls):   {total_neg:,}")
    print(f"Positive rate: {total_pos / len(combined):.3f}")


if __name__ == "__main__":
    main()
