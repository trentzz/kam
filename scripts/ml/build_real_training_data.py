#!/usr/bin/env python3
"""Extract 49-feature rows and labels from real Twist duplex TSV files.

Reads TSV outputs from the tsvs_disc_ml directory, labels each PASS row
against the truth VCF, applies the diff_pos -> VCF coordinate algorithm to
resolve genomic positions, and writes compressed feature CSVs for training
and testing.

Train split: ng 5 + ng 15.
Test split:  ng 30.

Outputs:
  output_dir/train_features.csv.gz
  output_dir/test_features.csv.gz
  output_dir/label_summary.csv

Usage:
    python3 scripts/ml/build_real_training_data.py \\
        --tsv-dir bigdata/experiments/03-ml-twist-duplex-v2/tsvs_disc_ml \\
        --truth-vcf docs/benchmarking/01-snvindel/scripts/truth_variants.vcf \\
        --output-dir bigdata/experiments/03-ml-twist-duplex-v2/

Run from the repository root.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Constants ──────────────────────────────────────────────────────────────────

# Variant class encoding — must match Rust ModelMeta.
VARIANT_CLASS_MAP: dict[str, int] = {
    "SNV": 0,
    "Insertion": 1,
    "Deletion": 2,
    "MNV": 3,
    "Complex": 4,
}

# Base encoding for trinucleotide context.
BASE_ENC: dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}

# Canonical substitution type encoding.
# Canonical classes pool complementary pairs.
SUBST_MAP: dict[tuple[str, str], int] = {
    ("C", "A"): 0,
    ("C", "G"): 1,
    ("C", "T"): 2,
    ("T", "A"): 3,
    ("T", "C"): 4,
    ("T", "G"): 5,
    ("G", "T"): 6,
    ("G", "C"): 7,
    ("G", "A"): 8,
    ("A", "T"): 9,
    ("A", "G"): 10,
    ("A", "C"): 11,
}

# The 49 feature names in order. The training/test CSVs use these as column
# headers; the ONNX model metadata records them for inference.
FEATURE_NAMES: list[str] = [
    # Original 33
    "vaf",
    "nref",
    "nalt",
    "ndupalt",
    "nsimalt",
    "sbp",
    "conf",
    "ref_len",
    "alt_len",
    "duplex_frac",
    "has_duplex",
    "ci_width",
    "alt_depth",
    "log_nalt",
    "log_nref",
    "log_alt_depth",
    "log_vaf",
    "vaf_times_conf",
    "vaf_times_nalt",
    "nalt_over_conf",
    "ci_width_rel",
    "snr",
    "conf_sq",
    "nalt_sq",
    "vaf_sq",
    "ref_alt_len_ratio",
    "indel_size",
    "duplex_enrichment",
    "simplex_only_frac",
    "conf_above_99",
    "conf_above_999",
    "sbp_above_05",
    "variant_class_enc",
    # Category B: 7 new TSV columns
    "n_simplex_fwd_alt",
    "n_simplex_rev_alt",
    "n_duplex_ref",
    "n_simplex_ref",
    "mean_alt_error_prob",
    "min_variant_specific_duplex",
    "mean_variant_specific_molecules",
    # Category C: 4 derived duplex/simplex features
    "strand_asymmetry_alt",
    "duplex_vaf",
    "simplex_vaf",
    "duplex_simplex_vaf_delta",
    # Category A: 7 sequence-context features
    "subst_type",
    "trinuc_context",
    "is_cpg",
    "gc_content_ref",
    "homopolymer_run",
    "dust_score",
    "repeat_fraction",
]


# ── Truth VCF parsing ──────────────────────────────────────────────────────────

def parse_truth_vcf(vcf_path: Path) -> set[tuple[str, int, str, str]]:
    """Parse a truth VCF and return a set of (chrom, pos, ref, alt) tuples.

    VCF positions are 1-based; they are kept as-is here because the
    diff_pos algorithm produces 0-based genomic positions relative to the
    target start. The final comparison uses VCF 1-based positions, so the
    truth set stores 1-based integer positions.

    Args:
        vcf_path: Path to the truth VCF file.

    Returns:
        Set of (chrom, genomic_pos_1based, ref, alt) tuples.

    Example:
        >>> s = parse_truth_vcf(Path("truth.vcf"))
        >>> ("chr1", 12345, "A", "T") in s
        False
    """
    truth: set[tuple[str, int, str, str]] = set()
    if not vcf_path.exists():
        print(f"[ERROR] Truth VCF not found: {vcf_path}", file=sys.stderr)
        return truth
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos_str, _id, ref, alt = (
                parts[0], parts[1], parts[2], parts[3], parts[4]
            )
            truth.add((chrom, int(pos_str), ref, alt))
    print(f"Loaded {len(truth):,} truth variants from {vcf_path.name}", flush=True)
    return truth


# ── Filename parsing ───────────────────────────────────────────────────────────

_FNAME_RE = re.compile(
    r"Sample_(\d+)ng_VAF_(\w+?)pc\.variants\.tsv$"
)


def parse_filename(fname: str) -> tuple[int, float] | None:
    """Extract (ng_condition, vaf_nominal) from a TSV filename.

    Handles VAF label encodings:
        0pc    -> 0.0
        0p001pc -> 0.001
        0p01pc  -> 0.01
        0p1pc   -> 0.1
        0p25pc  -> 0.25
        0p5pc   -> 0.5
        1pc     -> 1.0
        2pc     -> 2.0

    Args:
        fname: Basename of the TSV file, e.g. 'Sample_5ng_VAF_0p1pc.variants.tsv'.

    Returns:
        Tuple (ng_condition, vaf_nominal) or None if the filename does not match.

    Example:
        >>> parse_filename("Sample_5ng_VAF_0p1pc.variants.tsv")
        (5, 0.1)
        >>> parse_filename("Sample_30ng_VAF_0pc.variants.tsv")
        (30, 0.0)
    """
    m = _FNAME_RE.search(fname)
    if not m:
        return None
    ng = int(m.group(1))
    vaf_label = m.group(2)
    # Replace 'p' separator: '0p1' -> 0.1, '0p001' -> 0.001, '0' -> 0.0
    vaf_str = vaf_label.replace("p", ".")
    try:
        vaf = float(vaf_str)
    except ValueError:
        return None
    return ng, vaf


# ── diff_pos -> VCF coordinate algorithm ──────────────────────────────────────

def ref_alt_to_vcf_key(
    chrom: str,
    target_start: int,
    ref_seq: str,
    alt_seq: str,
) -> tuple[str, int, str, str] | None:
    """Convert a kam ref_seq/alt_seq pair to a VCF-style variant key.

    Uses the diff_pos algorithm from run_titration_batch.py lines 112-190.
    Produces a 1-based genomic position to match the truth VCF.

    Args:
        chrom: Chromosome name, e.g. 'chr1'.
        target_start: 0-based start of the target region.
        ref_seq: Full reference k-mer path from the TSV.
        alt_seq: Full alternate k-mer path from the TSV.

    Returns:
        Tuple (chrom, pos_1based, ref_allele, alt_allele) or None if the
        diff cannot be resolved.

    Example:
        >>> ref_alt_to_vcf_key("chr1", 100, "ACGT", "ATGT")
        ('chr1', 102, 'C', 'T')
    """
    diff_pos = next(
        (i for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a),
        None,
    )
    if diff_pos is None:
        return None

    # Trim common suffix.
    common_suffix = 0
    while (
        common_suffix < len(ref_seq) - diff_pos - 1
        and common_suffix < len(alt_seq) - diff_pos - 1
        and ref_seq[-(common_suffix + 1)] == alt_seq[-(common_suffix + 1)]
    ):
        common_suffix += 1

    if common_suffix > 0:
        ref_trimmed = ref_seq[diff_pos:-common_suffix]
        alt_trimmed = alt_seq[diff_pos:-common_suffix]
    else:
        ref_trimmed = ref_seq[diff_pos:]
        alt_trimmed = alt_seq[diff_pos:]

    if len(ref_seq) == len(alt_seq):
        # SNV or MNV path.
        ref_allele = ref_trimmed[0] if len(ref_trimmed) == 1 else ref_trimmed
        alt_allele = alt_trimmed[0] if len(alt_trimmed) == 1 else alt_trimmed
        genomic_pos = target_start + diff_pos
    else:
        # Indel path: further trim common prefix and suffix within trimmed seqs.
        inner_cs = 0
        while (
            inner_cs < len(ref_trimmed)
            and inner_cs < len(alt_trimmed)
            and ref_trimmed[-(inner_cs + 1)] == alt_trimmed[-(inner_cs + 1)]
        ):
            inner_cs += 1
        ref_min = ref_trimmed[:-inner_cs] if inner_cs > 0 else ref_trimmed
        alt_min = alt_trimmed[:-inner_cs] if inner_cs > 0 else alt_trimmed

        inner_cp = 0
        for rv, av in zip(ref_min, alt_min):
            if rv != av:
                break
            inner_cp += 1
        ref_min = ref_min[inner_cp:]
        alt_min = alt_min[inner_cp:]
        indel_start = diff_pos + inner_cp

        if len(ref_seq) > len(alt_seq):
            # Deletion.
            del_seq = ref_min
            anchor_pos = indel_start - 1
            while anchor_pos > 0 and del_seq and ref_seq[anchor_pos] == del_seq[-1]:
                del_seq = del_seq[-1:] + del_seq[:-1]
                anchor_pos -= 1
            if anchor_pos >= 0:
                anchor = ref_seq[anchor_pos]
                ref_allele = anchor + del_seq
                alt_allele = anchor
                genomic_pos = target_start + anchor_pos
            else:
                ref_allele = del_seq
                alt_allele = ""
                genomic_pos = target_start + indel_start
        else:
            # Insertion.
            ins_seq = alt_min
            anchor_pos = indel_start - 1
            while anchor_pos > 0 and ins_seq and ref_seq[anchor_pos] == ins_seq[-1]:
                ins_seq = ins_seq[-1:] + ins_seq[:-1]
                anchor_pos -= 1
            if anchor_pos >= 0:
                anchor = ref_seq[anchor_pos]
                ref_allele = anchor
                alt_allele = anchor + ins_seq
                genomic_pos = target_start + anchor_pos
            else:
                ref_allele = ""
                alt_allele = ins_seq
                genomic_pos = target_start + indel_start

    # target_start + diff_pos already gives the 1-based VCF position directly,
    # matching the convention used by run_titration_batch.py.
    return (chrom, genomic_pos, ref_allele, alt_allele)


# ── Sequence-context features (Category A) ────────────────────────────────────

def find_diff_pos(ref_seq: str, alt_seq: str) -> int | None:
    """Return the index of the first differing position between two sequences.

    Args:
        ref_seq: Reference sequence string.
        alt_seq: Alternate sequence string.

    Returns:
        Integer index of the first mismatch, or None if sequences are identical.

    Example:
        >>> find_diff_pos("ACGT", "ATGT")
        1
    """
    return next(
        (i for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a),
        None,
    )


def compute_subst_type(row: pd.Series) -> int:
    """Encode the substitution type for an SNV row.

    Uses the first differing position between ref_seq and alt_seq. Only
    applies to SNV rows (variant_type == 'SNV'); returns 12 for all others.

    Args:
        row: DataFrame row with 'variant_type', 'ref_seq', 'alt_seq' columns.

    Returns:
        Integer in 0-12 where 12 means non-SNV or unknown.

    Example:
        >>> import pandas as pd
        >>> compute_subst_type(pd.Series({"variant_type": "SNV", "ref_seq": "ACG", "alt_seq": "ATG"}))
        4
    """
    if row["variant_type"] != "SNV":
        return 12
    ref_seq = str(row["ref_seq"])
    alt_seq = str(row["alt_seq"])
    diff_pos = find_diff_pos(ref_seq, alt_seq)
    if diff_pos is None:
        return 12
    rb = ref_seq[diff_pos].upper()
    ab = alt_seq[diff_pos].upper()
    if rb not in "ACGT" or ab not in "ACGT":
        return 12
    return SUBST_MAP.get((rb, ab), 12)


def compute_trinuc_context(row: pd.Series) -> int:
    """Encode the trinucleotide context at the SNV position.

    Returns 64 for non-SNV rows, edge positions, or unknown bases.

    The encoding is: base_m1 * 16 + base_0 * 4 + base_p1, where each base
    is encoded A=0, C=1, G=2, T=3.

    Args:
        row: DataFrame row with 'variant_type', 'ref_seq', 'alt_seq' columns.

    Returns:
        Integer in 0-63, or 64 as a sentinel.

    Example:
        >>> import pandas as pd
        >>> compute_trinuc_context(pd.Series({"variant_type": "SNV", "ref_seq": "TACG", "alt_seq": "TATG"}))
        37
    """
    if row["variant_type"] != "SNV":
        return 64
    ref_seq = str(row["ref_seq"])
    alt_seq = str(row["alt_seq"])
    diff_pos = find_diff_pos(ref_seq, alt_seq)
    if diff_pos is None or diff_pos < 1 or diff_pos + 1 >= len(ref_seq):
        return 64
    b_m1 = ref_seq[diff_pos - 1].upper()
    b_0 = ref_seq[diff_pos].upper()
    b_p1 = ref_seq[diff_pos + 1].upper()
    if b_m1 not in BASE_ENC or b_0 not in BASE_ENC or b_p1 not in BASE_ENC:
        return 64
    return BASE_ENC[b_m1] * 16 + BASE_ENC[b_0] * 4 + BASE_ENC[b_p1]


def compute_is_cpg(row: pd.Series) -> int:
    """Return 1 if the SNV occurs at a CpG site, else 0.

    CpG criteria:
      - C>T substitution (subst_type == 2) with G immediately after the
        variant position in ref_seq.
      - G>A substitution (subst_type == 8) with C immediately before the
        variant position in ref_seq.

    Args:
        row: DataFrame row with 'variant_type', 'ref_seq', 'alt_seq',
             'subst_type' columns.

    Returns:
        1 or 0.

    Example:
        >>> import pandas as pd
        >>> compute_is_cpg(pd.Series({"variant_type": "SNV", "ref_seq": "ACG", "alt_seq": "ATG", "subst_type": 2}))
        1
    """
    if row["variant_type"] != "SNV":
        return 0
    subst = int(row["subst_type"])
    ref_seq = str(row["ref_seq"])
    alt_seq = str(row["alt_seq"])
    diff_pos = find_diff_pos(ref_seq, alt_seq)
    if diff_pos is None:
        return 0
    if subst == 2:
        # C>T: check the base immediately after diff_pos is G.
        if diff_pos + 1 < len(ref_seq) and ref_seq[diff_pos + 1].upper() == "G":
            return 1
    elif subst == 8:
        # G>A: check the base immediately before diff_pos is C.
        if diff_pos >= 1 and ref_seq[diff_pos - 1].upper() == "C":
            return 1
    return 0


def compute_gc_content(ref_seq: str) -> float:
    """Return the GC fraction of ref_seq.

    Args:
        ref_seq: Reference sequence string.

    Returns:
        Float in [0.0, 1.0]. Returns 0.0 for empty strings.

    Example:
        >>> compute_gc_content("ACGT")
        0.5
    """
    if not ref_seq:
        return 0.0
    gc = sum(1 for b in ref_seq.upper() if b in "GC")
    return gc / len(ref_seq)


def compute_homopolymer_run(row: pd.Series) -> int:
    """Return the longest run of identical bases adjacent to the variant position.

    Scans left from diff_pos-1 and right from diff_pos+1, counting matching
    bases. Returns the maximum of the two run lengths.

    Args:
        row: DataFrame row with 'ref_seq', 'alt_seq' columns.

    Returns:
        Integer run length (0 if diff_pos is at an edge or cannot be found).

    Example:
        >>> import pandas as pd
        >>> compute_homopolymer_run(pd.Series({"ref_seq": "AAACGT", "alt_seq": "AAATGT"}))
        3
    """
    ref_seq = str(row["ref_seq"])
    alt_seq = str(row["alt_seq"])
    diff_pos = find_diff_pos(ref_seq, alt_seq)
    if diff_pos is None:
        return 0
    # Left run: scan from diff_pos-1 backwards.
    left_run = 0
    if diff_pos >= 1:
        base = ref_seq[diff_pos - 1].upper()
        i = diff_pos - 1
        while i >= 0 and ref_seq[i].upper() == base:
            left_run += 1
            i -= 1
    # Right run: scan from diff_pos+1 forwards.
    right_run = 0
    if diff_pos + 1 < len(ref_seq):
        base = ref_seq[diff_pos + 1].upper()
        i = diff_pos + 1
        while i < len(ref_seq) and ref_seq[i].upper() == base:
            right_run += 1
            i += 1
    return max(left_run, right_run)


def compute_dust_score(ref_seq: str, window: int = 64) -> float:
    """Compute a DUST-like sequence complexity score.

    Higher values indicate more repetitive, low-complexity sequence.
    Uses trinucleotide frequency analysis over a sliding window.

    Args:
        ref_seq: Reference sequence string.
        window: Sliding window size in bases.

    Returns:
        Float >= 0.0. Values above 10 suggest borderline low complexity;
        above 30 indicates high repeat content.

    Example:
        >>> compute_dust_score("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") > 10
        True
        >>> compute_dust_score("ACGTGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG") < 10
        True
    """
    if not ref_seq or len(ref_seq) < 3:
        return 0.0
    seq = ref_seq.upper().encode("ascii", errors="replace")
    base_map = {ord(b"A"): 0, ord(b"C"): 1, ord(b"G"): 2, ord(b"T"): 3}
    max_score = 0.0
    win = min(window, len(seq))
    n_windows = max(1, len(seq) - win + 1)
    for start in range(n_windows):
        w = seq[start : start + win]
        counts: dict[tuple[int, int, int], int] = {}
        for i in range(len(w) - 2):
            tri = (
                base_map.get(w[i], 0),
                base_map.get(w[i + 1], 0),
                base_map.get(w[i + 2], 0),
            )
            counts[tri] = counts.get(tri, 0) + 1
        denom = max(1, len(w) - 2)
        score = sum(c * (c - 1) // 2 for c in counts.values()) / denom
        max_score = max(max_score, score)
    return max_score


def compute_repeat_fraction(ref_seq: str) -> float:
    """Compute fraction of bases in homopolymer (>=3) or dinucleotide (>=6) repeats.

    Args:
        ref_seq: Reference sequence string.

    Returns:
        Float in [0.0, 1.0]. Returns 0.0 for empty strings.

    Example:
        >>> compute_repeat_fraction("AAAAAACGT") > 0.5
        True
        >>> compute_repeat_fraction("ACGT")
        0.0
    """
    if not ref_seq:
        return 0.0
    seq = ref_seq.upper()
    n = len(seq)
    in_repeat = [False] * n

    # Homopolymer runs >= 3.
    i = 0
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if j - i >= 3:
            for k in range(i, j):
                in_repeat[k] = True
        i = j

    # Dinucleotide repeats >= 3 units (>= 6 bases).
    for start in range(n - 5):
        di = seq[start : start + 2]
        end = start + 2
        while end + 1 < n and seq[end] == di[0] and seq[end + 1] == di[1]:
            end += 2
        if end - start >= 6:
            for k in range(start, end):
                in_repeat[k] = True

    return sum(in_repeat) / n


# ── Feature computation ────────────────────────────────────────────────────────

def compute_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute all 51 features and return a new DataFrame with exactly those columns.

    The input DataFrame must have all required raw TSV columns. NaN values
    in numeric columns are filled with 0 before computation.

    Args:
        df: Raw TSV rows (PASS-filtered) with all expected columns.

    Returns:
        DataFrame with exactly the 51 columns listed in FEATURE_NAMES.

    Example:
        >>> import pandas as pd
        >>> row = pd.DataFrame([{"vaf": 0.01, "vaf_ci_low": 0.005, ...}])
        >>> feats = compute_features(row)
        >>> list(feats.columns) == FEATURE_NAMES
        True
    """
    out = pd.DataFrame(index=df.index)

    # ── Numeric inputs ────────────────────────────────────────────────────────
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

    ref_len = df["ref_seq"].fillna("").astype(str).str.len().astype(float)
    alt_len = df["alt_seq"].fillna("").astype(str).str.len().astype(float)
    ci_width = vaf_ci_high - vaf_ci_low
    alt_depth = nalt + nref

    # ── Original 33 features ──────────────────────────────────────────────────
    out["vaf"] = vaf
    out["nref"] = nref
    out["nalt"] = nalt
    out["ndupalt"] = ndupalt
    out["nsimalt"] = nsimalt
    out["sbp"] = sbp
    out["conf"] = conf
    out["ref_len"] = ref_len
    out["alt_len"] = alt_len
    out["duplex_frac"] = ndupalt / (nalt + 1e-9)
    out["has_duplex"] = (ndupalt > 0).astype(int)
    out["ci_width"] = ci_width
    out["alt_depth"] = alt_depth
    out["log_nalt"] = np.log1p(nalt)
    out["log_nref"] = np.log1p(nref)
    out["log_alt_depth"] = np.log1p(alt_depth)
    out["log_vaf"] = np.log(vaf + 1e-6)
    out["vaf_times_conf"] = vaf * conf
    out["vaf_times_nalt"] = vaf * nalt
    out["nalt_over_conf"] = nalt / (conf + 1e-9)
    out["ci_width_rel"] = ci_width / (vaf + 1e-9)
    out["snr"] = nalt / (nref + 1.0)
    out["conf_sq"] = conf ** 2
    out["nalt_sq"] = nalt ** 2
    out["vaf_sq"] = vaf ** 2
    out["ref_alt_len_ratio"] = ref_len / (alt_len + 1e-9)
    out["indel_size"] = (ref_len - alt_len).abs()
    # duplex_enrichment: how many duplex-alt molecules per expected duplex count.
    total_depth = nalt + nref
    out["duplex_enrichment"] = ndupalt / (vaf * total_depth + 1e-9)
    out["simplex_only_frac"] = nsimalt / (nalt + 1e-9)
    out["conf_above_99"] = (conf >= 0.99).astype(int)
    out["conf_above_999"] = (conf >= 0.999).astype(int)
    out["sbp_above_05"] = (sbp >= 0.05).astype(int)
    out["variant_class_enc"] = df["variant_type"].map(
        lambda v: VARIANT_CLASS_MAP.get(str(v), 0)
    )

    # ── Category B: 7 new TSV columns ─────────────────────────────────────────
    out["n_simplex_fwd_alt"] = n_sim_fwd
    out["n_simplex_rev_alt"] = n_sim_rev
    out["n_duplex_ref"] = n_dup_ref
    out["n_simplex_ref"] = n_sim_ref
    out["mean_alt_error_prob"] = mean_err
    out["min_variant_specific_duplex"] = min_vs_dup
    out["mean_variant_specific_molecules"] = mean_vs_mol

    # ── Category C: 4 derived duplex/simplex features ─────────────────────────
    out["strand_asymmetry_alt"] = (
        (n_sim_fwd - n_sim_rev) / (n_sim_fwd + n_sim_rev + 1e-9)
    )
    out["duplex_vaf"] = ndupalt / (ndupalt + n_dup_ref + 1e-9)
    out["simplex_vaf"] = nsimalt / (nsimalt + n_sim_ref + 1e-9)
    out["duplex_simplex_vaf_delta"] = out["duplex_vaf"] - out["simplex_vaf"]

    # ── Category A: 5 sequence-context features ───────────────────────────────
    # These require row-wise access to ref_seq/alt_seq strings.
    out["subst_type"] = df.apply(compute_subst_type, axis=1)
    out["trinuc_context"] = df.apply(compute_trinuc_context, axis=1)

    # is_cpg depends on subst_type, so attach it temporarily.
    df_tmp = df.copy()
    df_tmp["subst_type"] = out["subst_type"]
    out["is_cpg"] = df_tmp.apply(compute_is_cpg, axis=1)

    out["gc_content_ref"] = df["ref_seq"].fillna("").astype(str).map(compute_gc_content)
    out["homopolymer_run"] = df.apply(compute_homopolymer_run, axis=1)
    out["dust_score"] = df["ref_seq"].fillna("").astype(str).apply(compute_dust_score)
    out["repeat_fraction"] = (
        df["ref_seq"].fillna("").astype(str).apply(compute_repeat_fraction)
    )

    return out[FEATURE_NAMES]


# ── TSV processing ─────────────────────────────────────────────────────────────

_TARGET_RE = re.compile(r"^(chr\w+):(\d+)-(\d+)$")


def process_tsv(
    tsv_path: Path,
    ng_condition: int,
    vaf_nominal: float,
    truth_set: set[tuple[str, int, str, str]],
) -> pd.DataFrame:
    """Load a TSV, compute features, and assign labels for all PASS rows.

    Args:
        tsv_path: Path to the TSV file.
        ng_condition: Numeric ng value extracted from the filename (5, 15, or 30).
        vaf_nominal: Nominal VAF extracted from the filename.
        truth_set: Set of (chrom, pos_1based, ref, alt) tuples from the truth VCF.

    Returns:
        DataFrame with 49 feature columns plus label, sample_id, ng_condition,
        and vaf_nominal. Returns empty DataFrame if the file cannot be read.

    Example:
        >>> df = process_tsv(Path("sample.tsv"), 5, 0.001, set())
        >>> "label" in df.columns
        True
    """
    try:
        raw = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    except Exception as exc:
        print(f"[WARN] Could not read {tsv_path}: {exc}", file=sys.stderr)
        return pd.DataFrame()

    # Keep only PASS rows.
    if "filter" in raw.columns:
        raw = raw[raw["filter"] == "PASS"].copy()
    if raw.empty:
        return pd.DataFrame()

    # Assign labels.
    if vaf_nominal == 0.0:
        # Negative-only sample: every call is a false positive.
        raw["label"] = 0
    else:
        labels: list[int] = []
        for _, row in raw.iterrows():
            tid = str(row.get("target_id", ""))
            m = _TARGET_RE.match(tid)
            if not m:
                labels.append(0)
                continue
            chrom = m.group(1)
            target_start = int(m.group(2))
            ref_seq = str(row.get("ref_seq", ""))
            alt_seq = str(row.get("alt_seq", ""))
            vcf_key = ref_alt_to_vcf_key(chrom, target_start, ref_seq, alt_seq)
            if vcf_key is not None and vcf_key in truth_set:
                labels.append(1)
            else:
                labels.append(0)
        raw["label"] = labels

    feats = compute_features(raw)
    feats["label"] = raw["label"].values
    feats["sample_id"] = tsv_path.stem  # e.g. 'Sample_5ng_VAF_0p1pc.variants'
    feats["ng_condition"] = ng_condition
    feats["vaf_nominal"] = vaf_nominal
    return feats


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract 49-feature rows and labels from real Twist duplex TSV files.",
    )
    parser.add_argument(
        "--tsv-dir",
        type=Path,
        required=True,
        help="Directory containing Sample_*ng_VAF_*pc.variants.tsv files.",
    )
    parser.add_argument(
        "--truth-vcf",
        type=Path,
        required=True,
        help="Path to the truth VCF file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for train/test CSVs and label summary.",
    )
    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point."""
    args = parse_args()

    tsv_dir: Path = args.tsv_dir
    truth_vcf: Path = args.truth_vcf
    output_dir: Path = args.output_dir

    if not tsv_dir.exists():
        print(f"[ERROR] TSV directory not found: {tsv_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse truth VCF.
    truth_set = parse_truth_vcf(truth_vcf)

    # Discover and parse TSV files.
    tsv_files = sorted(tsv_dir.glob("Sample_*ng_VAF_*pc.variants.tsv"))
    if not tsv_files:
        print(f"[ERROR] No matching TSV files found in {tsv_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(tsv_files)} TSV files.", flush=True)

    train_frames: list[pd.DataFrame] = []
    test_frames: list[pd.DataFrame] = []
    summary_rows: list[dict] = []

    for tsv_path in tsv_files:
        parsed = parse_filename(tsv_path.name)
        if parsed is None:
            print(f"[WARN] Skipping unrecognised filename: {tsv_path.name}", file=sys.stderr)
            continue
        ng_condition, vaf_nominal = parsed

        df = process_tsv(tsv_path, ng_condition, vaf_nominal, truth_set)
        if df.empty:
            print(f"  [SKIP] {tsv_path.name}: no PASS rows", flush=True)
            continue

        n_total = len(df)
        n_tp = int((df["label"] == 1).sum())
        n_fp = int((df["label"] == 0).sum())
        print(
            f"  {tsv_path.name}: {n_total} rows, {n_tp} TP, {n_fp} FP",
            flush=True,
        )

        summary_rows.append({
            "sample_id": tsv_path.stem,
            "ng_condition": ng_condition,
            "vaf_nominal": vaf_nominal,
            "n_total": n_total,
            "n_tp": n_tp,
            "n_fp": n_fp,
        })

        # Train: ng 5 and ng 15. Test: ng 30.
        if ng_condition in (5, 15):
            train_frames.append(df)
        elif ng_condition == 30:
            test_frames.append(df)
        else:
            print(
                f"  [WARN] Unexpected ng={ng_condition}, skipping {tsv_path.name}",
                file=sys.stderr,
            )

    # Write outputs.
    col_order = FEATURE_NAMES + ["label", "sample_id", "ng_condition", "vaf_nominal"]

    if train_frames:
        train_df = pd.concat(train_frames, ignore_index=True)[col_order]
        train_path = output_dir / "train_features.csv.gz"
        train_df.to_csv(train_path, index=False, compression="gzip")
        n_pos = int((train_df["label"] == 1).sum())
        n_neg = int((train_df["label"] == 0).sum())
        print(
            f"\nTrain: {len(train_df):,} rows written to {train_path}",
            flush=True,
        )
        print(f"  Positives: {n_pos:,}  Negatives: {n_neg:,}", flush=True)
    else:
        print("[WARN] No train data collected.", file=sys.stderr)

    if test_frames:
        test_df = pd.concat(test_frames, ignore_index=True)[col_order]
        test_path = output_dir / "test_features.csv.gz"
        test_df.to_csv(test_path, index=False, compression="gzip")
        n_pos = int((test_df["label"] == 1).sum())
        n_neg = int((test_df["label"] == 0).sum())
        print(
            f"Test:  {len(test_df):,} rows written to {test_path}",
            flush=True,
        )
        print(f"  Positives: {n_pos:,}  Negatives: {n_neg:,}", flush=True)
    else:
        print("[WARN] No test data collected.", file=sys.stderr)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_dir / "label_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\nLabel summary written to {summary_path}", flush=True)

    # Per-sample breakdown.
    if not summary_df.empty:
        print("\nPer-sample breakdown:", flush=True)
        for _, row in summary_df.sort_values(["ng_condition", "vaf_nominal"]).iterrows():
            print(
                f"  {row['sample_id']:50s}  ng={row['ng_condition']:2d}  "
                f"vaf={row['vaf_nominal']:.4f}  "
                f"n={row['n_total']:5d}  tp={row['n_tp']:4d}  fp={row['n_fp']:5d}",
                flush=True,
            )


if __name__ == "__main__":
    main()
