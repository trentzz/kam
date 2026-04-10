"""Build ML training dataset for the twist-duplex experiment.

Reads simulated samples from:
  bigdata/experiments/01-ml-twist-duplex/simulations/train/<sample>/
  bigdata/experiments/01-ml-twist-duplex/simulations/test/<sample>/

Each sample directory must contain:
  calls_discovery.tsv        — discovery mode calls
  calls_tumour_informed.tsv  — tumour-informed mode calls
  params.json                — simulation parameters
  <sample>.truth.vcf.gz      — ground truth variants (gzipped VCF)

Outputs to bigdata/experiments/01-ml-twist-duplex/training_data.csv.
"""

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

TRAIN_ROOT = REPO_ROOT / "bigdata/experiments/01-ml-twist-duplex/simulations/train"
TEST_ROOT  = REPO_ROOT / "bigdata/experiments/01-ml-twist-duplex/simulations/test"

OUT_PATH = REPO_ROOT / "bigdata/experiments/01-ml-twist-duplex/training_data.csv"

# Rename TSV columns to the canonical CALL_COLS names.
COL_RENAME = {
    "target_id":       "chrom",
    "ref_seq":         "ref",
    "alt_seq":         "alt",
    "vaf_ci_low":      "vaf_lo",
    "vaf_ci_high":     "vaf_hi",
    "n_molecules_ref": "nref",
    "n_molecules_alt": "nalt",
    "n_duplex_alt":    "ndupalt",
    "n_simplex_alt":   "nsimalt",
    "strand_bias_p":   "sbp",
    "confidence":      "conf",
}

CALL_COLS = [
    "chrom", "variant_type", "ref", "alt", "filter", "vaf", "vaf_lo", "vaf_hi",
    "nref", "nalt", "ndupalt", "nsimalt", "sbp", "conf",
]

# Columns excluded from the output CSV. The ref/alt haplotype sequences can
# be several hundred characters each and make the CSV impractical to load for
# training. They are kept in CALL_COLS for labelling but dropped before output.
DROP_FROM_OUTPUT = {"ref", "alt", "chrom", "filter", "variant_class",
                    "dataset_type", "mode", "params_coverage",
                    "params_family_size_mean", "params_pcr_cycles",
                    "params_fragment_mean"}

TRUTH_KEY = ["chrom", "pos", "ref", "alt"]

PARAMS_COLS = [
    "params_coverage",
    "params_family_size_mean",
    "params_pcr_cycles",
    "params_fragment_mean",
]


def parse_ml3_chrom(chrom: str) -> tuple:
    """Parse a kam ML3 chrom descriptor to (chromosome, start, end).

    Format: 'chr1:49-250_SNV_1bp' → ('chr1', 49, 250).
    Plain chromosome: 'chr1' → ('chr1', None, None).
    """
    if ":" not in chrom:
        return chrom, None, None
    base, rest = chrom.split(":", 1)
    coords = rest.split("_")[0]  # '49-250'
    parts = coords.split("-")
    try:
        return base, int(parts[0]), int(parts[1])
    except (IndexError, ValueError):
        return base, None, None


def load_calls(path: Path) -> pd.DataFrame:
    """Load a calls TSV, applying column renames. Returns empty DataFrame on missing file."""
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns=COL_RENAME)
    for col in CALL_COLS:
        if col not in df.columns:
            df[col] = float("nan")
    return df[CALL_COLS].copy()


def load_truth_vcf(path: Path) -> pd.DataFrame:
    """Parse a gzipped VCF and return a DataFrame with chrom, pos, ref, alt.

    Multi-allelic ALT entries are reduced to the first allele.
    Returns an empty DataFrame if the file does not exist.
    """
    if not path.exists():
        return pd.DataFrame(columns=TRUTH_KEY)
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            pos   = int(parts[1])
            ref   = parts[3]
            alt   = parts[4].split(",")[0]  # first allele only
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt})
    if not rows:
        return pd.DataFrame(columns=TRUTH_KEY)
    return pd.DataFrame(rows)[TRUTH_KEY].drop_duplicates()


def load_params(sample_dir: Path) -> dict:
    """Load params.json and return flattened numeric fields plus split."""
    params_path = sample_dir / "params.json"
    if not params_path.exists():
        return {}
    with open(params_path) as f:
        raw = json.load(f)
    return {
        "params_coverage":         raw.get("coverage", float("nan")),
        "params_family_size_mean": raw.get("family_size_mean", float("nan")),
        "params_pcr_cycles":       raw.get("pcr_cycles", float("nan")),
        "params_fragment_mean":    raw.get("fragment_mean", float("nan")),
        "split":                   raw.get("split", "train"),
        "vtype":                   raw.get("vtype", ""),
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
    vt_str_map = {
        "SNV": "SNV", "MNV": "MNV",
        "Deletion": "Deletion", "Insertion": "Insertion",
        "LargeDeletion": "LargeDeletion", "TandemDuplication": "TandemDuplication",
        "Inversion": "Inversion", "InvDel": "InvDel",
        "NovelInsertion": "NovelInsertion", "Complex": "Complex", "Fusion": "Fusion",
    }

    if "variant_type" in df.columns and df["variant_type"].notna().any():
        df["variant_class"] = df["variant_type"].map(vt_str_map).fillna("SNV")
    else:
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
    tsv_name = "calls_discovery.tsv" if mode == "discovery" else "calls_tumour_informed.tsv"
    calls = load_calls(sample_dir / tsv_name)
    if calls.empty:
        return pd.DataFrame()

    # Build truth list as (chrom, pos, ref, alt) for exact position matching.
    truth_full = [
        (str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"]))
        for _, row in truth_df.iterrows()
        if not pd.isna(row.get("pos", float("nan")))
    ]

    is_invdel = params.get("vtype", "") == "invdel"

    def label_call(row: pd.Series) -> int:
        if is_invdel:
            # All calls in an invdel sample are positive if truth has variants.
            return 1 if truth_full else 0

        chrom_base, start, end = parse_ml3_chrom(str(row["chrom"]))
        if start is None:
            # Plain chrom format (not expected in twist-duplex data).
            return 1 if any(t == chrom_base for t, _, _, _ in truth_full) else 0

        # Check if any truth variant falls in this window.
        truth_in_window = [
            (t_chrom, t_pos)
            for t_chrom, t_pos, _, _ in truth_full
            if t_chrom == chrom_base and start <= t_pos <= end
        ]
        if not truth_in_window:
            return 0

        # Find the genomic position of the first sequence difference.
        # Calls encode full window haplotypes; the first diff gives the variant site.
        # Window start is 0-based; VCF positions are 1-based.
        ref_str = str(row["ref"])
        alt_str = str(row["alt"])
        first_diff_genomic = None
        for i, (r, a) in enumerate(zip(ref_str, alt_str)):
            if r != a:
                first_diff_genomic = start + i + 1  # convert to 1-based VCF coords
                break

        if first_diff_genomic is None:
            return 0

        # Label as positive if the variant site maps to a truth position.
        # VCF anchor-base convention: deletions/insertions record the anchor position,
        # but the sequence diff appears at anchor+1. Check both positions.
        for t_chrom, t_pos in truth_in_window:
            if first_diff_genomic == t_pos or first_diff_genomic == t_pos + 1:
                return 1
        return 0

    calls["label"] = calls.apply(label_call, axis=1)
    calls = derive_features(calls)
    calls["sample_id"]    = sample_id
    calls["dataset_type"] = dataset_type
    calls["mode"]         = mode
    calls["split"]        = split
    for col in PARAMS_COLS:
        calls[col] = params.get(col, float("nan"))
    # Drop large string columns not needed for training.
    drop = [c for c in DROP_FROM_OUTPUT if c in calls.columns]
    return calls.drop(columns=drop)


def collect_samples(
    sample_root: Path,
    dataset_type: str,
    split: str,
    out_path: Path,
    first_chunk: bool,
) -> tuple[int, int]:
    """Process all sample directories under sample_root, writing rows to out_path.

    Writes in chunks to keep memory usage flat. Returns (n_rows_written, n_skipped).
    """
    if not sample_root.exists():
        return 0, 0

    n_rows = 0
    n_skip = 0
    chunk: list[pd.DataFrame] = []
    chunk_rows = 0
    chunk_limit = 50_000  # flush every ~50k rows

    sample_dirs = sorted(p for p in sample_root.iterdir() if p.is_dir())
    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        truth_path = sample_dir / f"{sample_id}.truth.vcf.gz"
        truth_df = load_truth_vcf(truth_path)
        if truth_df.empty:
            n_skip += 1
            continue
        params = load_params(sample_dir)
        effective_split = params.get("split", split)
        for mode in ("discovery", "tumour_informed"):
            df = process_mode(
                sample_dir, sample_id, dataset_type, mode,
                truth_df, params, effective_split,
            )
            if not df.empty:
                chunk.append(df)
                chunk_rows += len(df)

        if chunk_rows >= chunk_limit:
            combined = pd.concat(chunk, ignore_index=True)
            write_header = first_chunk
            combined.to_csv(out_path, mode="a", header=write_header, index=False)
            n_rows += len(combined)
            first_chunk = False
            chunk = []
            chunk_rows = 0

    # Flush remaining.
    if chunk:
        combined = pd.concat(chunk, ignore_index=True)
        combined.to_csv(out_path, mode="a", header=first_chunk, index=False)
        n_rows += len(combined)

    return n_rows, n_skip


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    # Truncate any existing file.
    OUT_PATH.write_text("")

    n_train, skip_train = collect_samples(TRAIN_ROOT, "twist_duplex", "train", OUT_PATH, first_chunk=True)
    n_test,  skip_test  = collect_samples(TEST_ROOT,  "twist_duplex", "test",  OUT_PATH, first_chunk=(n_train == 0))

    total = n_train + n_test
    if total == 0:
        print("ERROR: no data collected", file=sys.stderr)
        sys.exit(1)

    print(f"Written {total:,} rows to {OUT_PATH}")

    # Read just the split/label columns to print the distribution.
    summary = pd.read_csv(OUT_PATH, usecols=["split", "label"])
    print(f"Columns: {len(pd.read_csv(OUT_PATH, nrows=0).columns)}")
    print("\nSplit distribution:")
    split_stats = (
        summary.groupby(["split", "label"])
        .size()
        .rename("count")
        .reset_index()
    )
    print(split_stats.to_string(index=False))


if __name__ == "__main__":
    main()
