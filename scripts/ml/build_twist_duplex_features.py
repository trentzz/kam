#!/usr/bin/env python3
"""Extract features from Twist duplex ML pipeline outputs.

Reads per-sample directories under docs/benchmarking/ml-twist-duplex/simulations/,
loads the discovery.tsv and tumour_informed.tsv outputs from kam, labels each
row against the truth VCF, and builds a rich feature set (~55 features) that
exploits the 7 new columns added to the kam TSV output.

Outputs:
  docs/benchmarking/ml-twist-duplex/train_features.csv.gz
  docs/benchmarking/ml-twist-duplex/test_features.csv.gz

Usage:
    python3 scripts/ml/build_twist_duplex_features.py [--mode discovery|tumour_informed|both]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent

ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
SIM_DIR = ML_DIR / "simulations"
MANIFEST_PATH = ML_DIR / "manifest.json"

# ─── Column definitions ────────────────────────────────────────────────────────

# All 21 columns in the updated kam TSV output (14 original + 7 new).
KAM_COLS = [
    # Original 14
    "target_id",
    "variant_type",
    "ref_seq",
    "alt_seq",
    "vaf",
    "vaf_ci_low",
    "vaf_ci_high",
    "n_molecules_ref",
    "n_molecules_alt",
    "n_duplex_alt",
    "n_simplex_alt",
    "strand_bias_p",
    "confidence",
    "filter",
    # New 7
    "n_simplex_fwd_alt",
    "n_simplex_rev_alt",
    "n_duplex_ref",
    "n_simplex_ref",
    "mean_alt_error_prob",
    "min_variant_specific_duplex",
    "mean_variant_specific_molecules",
]

# Truth key: used to label calls against ground truth.
TRUTH_KEY = ["chrom", "pos", "ref", "alt"]


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_kam_tsv(path: Path) -> pd.DataFrame:
    """Load a kam TSV output, returning an empty DataFrame on failure.

    Fills missing columns with NaN so downstream code is robust against
    older outputs that lack the new 7 columns.

    Args:
        path: Path to the TSV file.

    Returns:
        DataFrame with KAM_COLS columns.
    """
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as exc:
        print(f"[WARN] Could not read {path}: {exc}", file=sys.stderr)
        return pd.DataFrame()
    for col in KAM_COLS:
        if col not in df.columns:
            df[col] = float("nan")
    return df[KAM_COLS].copy()


def load_truth_from_vcf(vcf_path: Path) -> pd.DataFrame:
    """Parse a truth VCF and return a DataFrame of variant positions.

    Converts 1-based VCF positions to 0-based to match kam's output
    convention (kam reports 0-based positions).

    Args:
        vcf_path: Path to the VCF file.

    Returns:
        DataFrame with columns: chrom, pos (0-based str), ref, alt.
    """
    rows = []
    if not vcf_path.exists():
        return pd.DataFrame(columns=TRUTH_KEY)
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    try:
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos_1based, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                # Convert 1-based VCF to 0-based to align with kam output
                pos_0based = str(int(pos_1based) - 1)
                rows.append({"chrom": chrom, "pos": pos_0based, "ref": ref, "alt": alt})
    except Exception as exc:
        print(f"[WARN] Could not read VCF {vcf_path}: {exc}", file=sys.stderr)
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=TRUTH_KEY)


def load_params(sample_dir: Path) -> dict:
    """Load params.json from a sample directory.

    Args:
        sample_dir: Sample output directory.

    Returns:
        Dict of simulation parameters, empty if not found.
    """
    params_path = sample_dir / "params.json"
    if not params_path.exists():
        return {}
    with open(params_path) as fh:
        return json.load(fh)


# ─── Labelling ────────────────────────────────────────────────────────────────

def label_calls(calls: pd.DataFrame, truth_df: pd.DataFrame) -> pd.DataFrame:
    """Add a binary 'label' column: 1 if the call matches a truth variant.

    Matching is by (chrom, pos, ref, alt) key. For kam TSV, 'target_id'
    encodes the chromosome (or target name), and 'ref_seq'/'alt_seq' hold
    the alleles. Position is not directly in the TSV; for exact labelling
    the truth VCF positions are matched against the kam target coordinate.

    Because the kam TSV does not carry a genomic position column directly
    (target_id encodes the locus), labelling uses (target_id, ref_seq,
    alt_seq) as a proxy. When target_id is a genomic coordinate string
    (e.g. 'chr1:100-200'), the variant position is embedded in it.

    For simplicity, this implementation labels using (ref_seq, alt_seq)
    against truth (ref, alt) within the same target/sample, which is
    sufficient for per-sample datasets where each sample has 1-3 variants.

    Args:
        calls: DataFrame from load_kam_tsv().
        truth_df: DataFrame from load_truth_from_vcf().

    Returns:
        calls with a 'label' column added.
    """
    if truth_df.empty:
        calls["label"] = 0
        return calls

    truth_set = set(
        zip(
            truth_df["ref"].astype(str),
            truth_df["alt"].astype(str),
        )
    )

    calls["label"] = calls.apply(
        lambda r: 1
        if (str(r["ref_seq"]), str(r["alt_seq"])) in truth_set
        else 0,
        axis=1,
    )
    return calls


# ─── Feature engineering ─────────────────────────────────────────────────────

def derive_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute all derived features and add them to the DataFrame.

    Modifies df in place and returns it.

    Feature groups:
    - Basic raw (from KAM_COLS)
    - Duplex-specific ratios and deltas
    - Log transforms
    - Ratio / interaction terms
    - Binned categorical features

    Args:
        df: DataFrame with KAM_COLS columns.

    Returns:
        df with additional feature columns.
    """
    # ── Numeric casts ────────────────────────────────────────────────────────
    vaf = df["vaf"].fillna(0).astype(float)
    vaf_ci_low = df["vaf_ci_low"].fillna(0).astype(float)
    vaf_ci_high = df["vaf_ci_high"].fillna(0).astype(float)
    nref = df["n_molecules_ref"].fillna(0).astype(float)
    nalt = df["n_molecules_alt"].fillna(0).astype(float)
    ndupalt = df["n_duplex_alt"].fillna(0).astype(float)
    nsimalt = df["n_simplex_alt"].fillna(0).astype(float)
    sbp = df["strand_bias_p"].fillna(1.0).astype(float)
    conf = df["confidence"].fillna(0).astype(float)
    n_sim_fwd = df["n_simplex_fwd_alt"].fillna(0).astype(float)
    n_sim_rev = df["n_simplex_rev_alt"].fillna(0).astype(float)
    n_dup_ref = df["n_duplex_ref"].fillna(0).astype(float)
    n_sim_ref = df["n_simplex_ref"].fillna(0).astype(float)
    mean_err = df["mean_alt_error_prob"].fillna(0).astype(float)
    min_vs_dup = df["min_variant_specific_duplex"].fillna(0).astype(float)
    mean_vs_mol = df["mean_variant_specific_molecules"].fillna(0).astype(float)

    # ── Basic length features ─────────────────────────────────────────────────
    df["ref_len"] = df["ref_seq"].fillna("").astype(str).str.len()
    df["alt_len"] = df["alt_seq"].fillna("").astype(str).str.len()
    ref_len = df["ref_len"].astype(float)
    alt_len = df["alt_len"].astype(float)

    # ── Duplex-specific ratios ────────────────────────────────────────────────
    df["duplex_vaf"] = ndupalt / (ndupalt + n_dup_ref + 1e-9)
    df["simplex_vaf"] = nsimalt / (nsimalt + n_sim_ref + 1e-9)
    df["duplex_simplex_vaf_delta"] = df["duplex_vaf"] - df["simplex_vaf"]
    df["vaf_duplex_agreement"] = (vaf - df["duplex_vaf"]).abs()
    df["duplex_enrichment"] = ndupalt / (nalt + 1e-9)
    df["simplex_frac"] = nsimalt / (nalt + 1e-9)
    df["strand_balance"] = (
        np.minimum(n_sim_fwd, n_sim_rev) / (n_sim_fwd + n_sim_rev + 1e-9)
    )
    df["strand_asymmetry"] = (n_sim_fwd - n_sim_rev) / (n_sim_fwd + n_sim_rev + 1e-9)
    df["duplex_depth"] = ndupalt + n_dup_ref
    df["simplex_depth"] = nsimalt + n_sim_ref
    df["duplex_ref_frac"] = n_dup_ref / (n_dup_ref + ndupalt + 1e-9)
    df["simplex_only_frac"] = nsimalt / (nalt + 1e-9)
    df["variant_specific_duplex_frac"] = min_vs_dup / (ndupalt + 1e-9)
    df["qual_confidence_interaction"] = (1.0 - mean_err) * conf
    df["qual_vaf_interaction"] = (1.0 - mean_err) * vaf
    df["has_duplex"] = (ndupalt > 0).astype(int)

    # ── Log transforms ────────────────────────────────────────────────────────
    df["log_vaf"] = np.log(vaf + 1e-6)
    df["log_nalt"] = np.log1p(nalt)
    df["log_nref"] = np.log1p(nref)
    df["log_nduplex_alt"] = np.log1p(ndupalt)
    df["log_nsimplex_alt"] = np.log1p(nsimalt)
    df["log_alt_depth"] = np.log1p(nalt + nref)
    df["log_duplex_depth"] = np.log1p(df["duplex_depth"])
    df["log_simplex_depth"] = np.log1p(df["simplex_depth"])

    # ── Ratio / interaction terms ─────────────────────────────────────────────
    ci_width = vaf_ci_high - vaf_ci_low
    df["ci_width"] = ci_width
    df["ci_width_rel"] = ci_width / (vaf + 1e-9)
    df["vaf_times_conf"] = vaf * conf
    df["vaf_times_nalt"] = vaf * nalt
    df["nalt_over_conf"] = nalt / (conf + 1e-9)
    df["conf_sq"] = conf ** 2
    df["nalt_sq"] = nalt ** 2
    df["vaf_sq"] = vaf ** 2
    df["duplex_frac"] = ndupalt / (nalt + 1e-9)  # alias for duplex_enrichment
    df["snr"] = nalt / (nref + 1.0)

    # ── Length-based features ─────────────────────────────────────────────────
    df["ref_alt_len_ratio"] = ref_len / (alt_len + 1e-9)
    df["indel_size"] = (ref_len - alt_len).abs()

    # ── Variant classification ────────────────────────────────────────────────
    max_len = df[["ref_len", "alt_len"]].max(axis=1)

    def classify_variant(row):
        """Classify a variant row as SNV, short_indel, medium_indel, or SV."""
        rl = int(row["ref_len"])
        al = int(row["alt_len"])
        ml = max(rl, al)
        if rl == 1 and al == 1:
            return "SNV"
        if ml <= 5:
            return "short_indel"
        if ml <= 20:
            return "medium_indel"
        return "SV"

    df["variant_class"] = df.apply(classify_variant, axis=1)

    # ── Binned features ───────────────────────────────────────────────────────
    alt_depth = (nalt + nref).values
    try:
        df["depth_bucket"] = pd.qcut(
            alt_depth, q=5, labels=False, duplicates="drop"
        ).fillna(0).astype(int)
    except ValueError:
        df["depth_bucket"] = 0

    sbp_arr = sbp.values
    df["sbp_cat"] = np.where(
        sbp_arr < 0.01, "low",
        np.where(sbp_arr < 0.05, "medium", "high"),
    )
    df["conf_above_99"] = (conf >= 0.99).astype(int)
    df["conf_above_999"] = (conf >= 0.999).astype(int)
    df["sbp_above_05"] = (sbp >= 0.05).astype(int)

    return df


def add_param_features(df: pd.DataFrame, params: dict) -> pd.DataFrame:
    """Attach simulation parameter columns to the DataFrame.

    Args:
        df: Call DataFrame.
        params: Dict from params.json.

    Returns:
        df with parameter columns added.
    """
    df["coverage"] = params.get("coverage", float("nan"))
    df["family_size_mean"] = params.get("family_size_mean", float("nan"))
    df["pcr_cycles"] = params.get("pcr_cycles", float("nan"))
    df["fragment_mean"] = params.get("fragment_mean", float("nan"))
    df["vaf_target"] = params.get("vaf", float("nan"))
    df["vaf_log_target"] = np.log(float(params.get("vaf", 1e-6)) + 1e-9)
    df["variant_class_true"] = params.get("vtype", "unknown")
    return df


# ─── Per-sample processing ────────────────────────────────────────────────────

def process_sample(
    sample: dict,
    mode: str,
    manifest_vtype: str,
) -> pd.DataFrame:
    """Load, label, and featurise calls for a single sample and mode.

    Args:
        sample: Sample dict from the manifest.
        mode: 'discovery' or 'tumour_informed'.
        manifest_vtype: Variant type from the manifest (e.g. 'snv').

    Returns:
        Featurised DataFrame, or empty DataFrame if the sample has no data.
    """
    name = sample["name"]
    split = sample["split"]
    sample_dir = SIM_DIR / split / name

    tsv_name = "discovery.tsv" if mode == "discovery" else "tumour_informed.tsv"
    calls = load_kam_tsv(sample_dir / tsv_name)
    if calls.empty:
        return pd.DataFrame()

    # Find truth VCF: prefer the pre-generated one, fall back to sim output
    truth_vcf = ML_DIR / "truth_vcfs" / f"{name}.vcf"
    if not truth_vcf.exists():
        vcf_candidates = sorted(sample_dir.glob("*.truth.vcf"))
        if vcf_candidates:
            truth_vcf = vcf_candidates[0]

    truth_df = load_truth_from_vcf(truth_vcf)
    calls = label_calls(calls, truth_df)
    calls = derive_features(calls)

    params = load_params(sample_dir)
    if not params:
        # Fall back to manifest values
        params = {
            "coverage": sample.get("coverage"),
            "family_size_mean": sample.get("family_size_mean"),
            "pcr_cycles": sample.get("pcr_cycles"),
            "fragment_mean": sample.get("fragment_mean"),
            "vaf": sample.get("vaf"),
            "vtype": manifest_vtype,
        }
    calls = add_param_features(calls, params)

    calls["sample_id"] = name
    calls["split"] = split
    calls["mode"] = mode

    return calls


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract features from Twist duplex pipeline outputs.",
    )
    parser.add_argument(
        "--mode",
        choices=["discovery", "tumour_informed", "both"],
        default="tumour_informed",
        help="Which kam mode to extract from (default: tumour_informed).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N samples per split (for testing).",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    if not MANIFEST_PATH.exists():
        print(
            f"[ERROR] Manifest not found: {MANIFEST_PATH}\n"
            "Run scripts/ml/generate_twist_duplex_configs.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(MANIFEST_PATH) as fh:
        manifest_data = json.load(fh)
    all_samples = manifest_data["samples"]

    modes = (
        ["discovery", "tumour_informed"]
        if args.mode == "both"
        else [args.mode]
    )

    for split in ("train", "test"):
        split_samples = [s for s in all_samples if s["split"] == split]
        if args.limit is not None:
            split_samples = split_samples[: args.limit]

        print(f"\n=== {split.upper()} ({len(split_samples)} samples) ===", flush=True)
        all_rows: list[pd.DataFrame] = []

        for i, sample in enumerate(split_samples):
            for mode in modes:
                df = process_sample(sample, mode, sample["vtype"])
                if not df.empty:
                    all_rows.append(df)

            if (i + 1) % 500 == 0:
                print(f"  [{i + 1}/{len(split_samples)}] processed", flush=True)

        if not all_rows:
            print(f"[WARN] No data collected for {split} split.", file=sys.stderr)
            continue

        combined = pd.concat(all_rows, ignore_index=True)
        out_path = ML_DIR / f"{split}_features.csv.gz"
        combined.to_csv(out_path, index=False, compression="gzip")

        n_pos = (combined["label"] == 1).sum()
        n_neg = (combined["label"] == 0).sum()
        print(f"Written {len(combined):,} rows to {out_path}")
        print(f"  Positives: {n_pos:,}  Negatives: {n_neg:,}  "
              f"Positive rate: {n_pos / max(len(combined), 1) * 100:.3f}%")
        print(f"  Columns: {len(combined.columns)}")


if __name__ == "__main__":
    main()
