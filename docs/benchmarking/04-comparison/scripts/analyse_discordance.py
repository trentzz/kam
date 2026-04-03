"""
analyse_discordance.py

Analyse discordant variant calls between alignment-based and kam methods.

Reads concordance.csv (produced by build_concordance.py) and classifies each
discordant row by the most likely reason for the discordance.

Discordance directions
----------------------
  alignment_only — alignment detected the variant; kam did not.
  kam_only       — kam detected the variant; alignment did not.

Reason categories
-----------------
  low_depth      — alignment_depth < LOW_DEPTH_THRESHOLD (or no depth recorded).
  low_vaf        — vaf_level is at or below the near-detection-limit threshold.
  variant_type   — discordance clusters on a specific variant type (INDEL / SV).
  input_amount   — discordance clusters at a specific input_ng level.
  unknown        — does not fit the above categories.

The categories are applied in priority order: a row can only belong to one
category (the first that matches).

Outputs
-------
  discordance_details.csv   — one row per discordant variant with reason.
  discordance_summary.csv   — counts by (category, direction, variant_type).
  Stdout summary.

Usage
-----
  python analyse_discordance.py [--concordance PATH] [--out-dir PATH]
                                [--depth-threshold INT] [--low-vaf-labels LABEL...]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
COMPARISON_DIR = SCRIPT_DIR.parent

DEFAULT_CONCORDANCE = COMPARISON_DIR / "concordance.csv"
DEFAULT_OUT_DIR = COMPARISON_DIR

# Molecule depth below this value is classified as low_depth.
DEFAULT_DEPTH_THRESHOLD = 10

# VAF labels considered near the detection limit.
DEFAULT_LOW_VAF_LABELS = {"0p001pc", "0p01pc", "0p1pc", "0p25pc"}

# Variant types that commonly cause type-specific discordance.
TYPE_SPECIFIC_TYPES = {"INDEL", "INS", "DEL", "SV", "INVDEL"}

# Input amounts (ng) considered low-input.
LOW_INPUT_NG = {5}


# ---------------------------------------------------------------------------
# Reason classification
# ---------------------------------------------------------------------------

def classify_reason(
    row: dict,
    depth_threshold: int,
    low_vaf_labels: set[str],
) -> str:
    """Return the most likely reason category for a discordant variant.

    Priority order:
      1. low_depth  — alignment_depth is present and below threshold, or absent
                      for alignment_only variants (depth not recorded implies
                      the variant was not covered adequately).
      2. low_vaf    — vaf_level is in the near-detection-limit set.
      3. variant_type — variant_type is an INDEL or structural type.
      4. input_amount — input_ng is in the low-input set.
      5. unknown    — none of the above apply.
    """
    direction = row["concordance"]

    # 1. Low depth.
    depth_str = row.get("alignment_depth", "").strip()
    if direction == "alignment_only":
        # alignment detected it; check whether depth was very low (possibly a
        # spurious alignment call at borderline depth).
        if depth_str:
            try:
                depth = int(depth_str)
                if depth < depth_threshold:
                    return "low_depth"
            except ValueError:
                pass
        else:
            # No depth recorded alongside a detection — treat as low depth.
            return "low_depth"
    else:
        # kam_only: alignment did not detect it; missing/low depth explains why.
        if not depth_str:
            return "low_depth"
        try:
            depth = int(depth_str)
            if depth < depth_threshold:
                return "low_depth"
        except ValueError:
            return "low_depth"

    # 2. Low VAF.
    if row.get("vaf_level", "") in low_vaf_labels:
        return "low_vaf"

    # 3. Variant type.
    if row.get("variant_type", "") in TYPE_SPECIFIC_TYPES:
        return "variant_type"

    # 4. Input amount.
    try:
        input_ng = int(row.get("input_ng", "0"))
        if input_ng in LOW_INPUT_NG:
            return "input_amount"
    except ValueError:
        pass

    return "unknown"


# ---------------------------------------------------------------------------
# Load concordance
# ---------------------------------------------------------------------------

def load_concordance(path: Path) -> list[dict]:
    """Load concordance.csv and return all rows as a list of dicts."""
    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Build discordance details
# ---------------------------------------------------------------------------

def build_discordance_details(
    concordance_rows: list[dict],
    depth_threshold: int,
    low_vaf_labels: set[str],
) -> list[dict]:
    """Extract discordant rows and annotate with reason and direction-specific fields.

    Returns a list of dicts with the following columns:
      sample_id, input_ng, vaf_level, chromosome, position, ref, alt,
      variant_type, concordance (direction), alignment_vaf, alignment_depth,
      kam_vaf, kam_confidence, reason
    """
    details = []

    for row in concordance_rows:
        direction = row.get("concordance", "")
        if direction not in ("alignment_only", "kam_only"):
            continue

        reason = classify_reason(row, depth_threshold, low_vaf_labels)

        details.append({
            "sample_id": row["sample_id"],
            "input_ng": row["input_ng"],
            "vaf_level": row["vaf_level"],
            "chromosome": row["chromosome"],
            "position": row["position"],
            "ref": row["ref"],
            "alt": row["alt"],
            "variant_type": row["variant_type"],
            "direction": direction,
            "alignment_vaf": row.get("alignment_vaf", ""),
            "alignment_depth": row.get("alignment_depth", ""),
            "kam_vaf": row.get("kam_vaf", ""),
            "kam_confidence": row.get("kam_confidence", ""),
            "reason": reason,
        })

    # Deterministic sort.
    details.sort(key=lambda r: (r["direction"], r["reason"], r["vaf_level"], r["chromosome"], int(r["position"])))

    return details


# ---------------------------------------------------------------------------
# Build discordance summary
# ---------------------------------------------------------------------------

def build_discordance_summary(details: list[dict]) -> list[dict]:
    """Count discordant variants by (reason, direction, variant_type).

    Returns rows sorted by direction, then reason, then variant_type.
    """
    counts: dict[tuple[str, str, str], int] = defaultdict(int)

    for row in details:
        key = (row["reason"], row["direction"], row["variant_type"])
        counts[key] += 1

    summary = []
    for (reason, direction, vtype), count in sorted(counts.items()):
        summary.append({
            "reason": reason,
            "direction": direction,
            "variant_type": vtype,
            "count": count,
        })

    return summary


# ---------------------------------------------------------------------------
# Print stdout summary
# ---------------------------------------------------------------------------

def print_summary(details: list[dict], depth_threshold: int) -> None:
    """Print a human-readable summary of discordance to stdout."""
    total = len(details)
    if total == 0:
        print("No discordant variants found.")
        return

    align_only = sum(1 for r in details if r["direction"] == "alignment_only")
    kam_only = sum(1 for r in details if r["direction"] == "kam_only")

    reason_counts: dict[str, int] = defaultdict(int)
    for r in details:
        reason_counts[r["reason"]] += 1

    vtype_counts: dict[str, int] = defaultdict(int)
    for r in details:
        vtype_counts[r["variant_type"]] += 1

    print(f"Discordance summary ({total} discordant variants)")
    print(f"  depth threshold used : {depth_threshold} molecules")
    print()
    print("  Direction:")
    print(f"    alignment_only  : {align_only:>5}  ({100 * align_only / total:.1f}%)")
    print(f"    kam_only        : {kam_only:>5}  ({100 * kam_only / total:.1f}%)")
    print()
    print("  Reason breakdown:")
    for reason in ("low_depth", "low_vaf", "variant_type", "input_amount", "unknown"):
        c = reason_counts.get(reason, 0)
        print(f"    {reason:<16}: {c:>5}  ({100 * c / total:.1f}%)")
    print()
    print("  Variant type breakdown:")
    for vtype, c in sorted(vtype_counts.items()):
        print(f"    {vtype:<16}: {c:>5}  ({100 * c / total:.1f}%)")


# ---------------------------------------------------------------------------
# Field lists
# ---------------------------------------------------------------------------

DETAILS_FIELDS = [
    "sample_id",
    "input_ng",
    "vaf_level",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant_type",
    "direction",
    "alignment_vaf",
    "alignment_depth",
    "kam_vaf",
    "kam_confidence",
    "reason",
]

SUMMARY_FIELDS = [
    "reason",
    "direction",
    "variant_type",
    "count",
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Analyse discordant variants from concordance.csv."
    )
    parser.add_argument(
        "--concordance",
        type=Path,
        default=DEFAULT_CONCORDANCE,
        help="Path to concordance.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for discordance CSV files (default: %(default)s)",
    )
    parser.add_argument(
        "--depth-threshold",
        type=int,
        default=DEFAULT_DEPTH_THRESHOLD,
        help="Molecule depth below which a variant is classified as low_depth "
             "(default: %(default)s)",
    )
    parser.add_argument(
        "--low-vaf-labels",
        nargs="+",
        default=sorted(DEFAULT_LOW_VAF_LABELS),
        help="VAF level labels considered near the detection limit "
             "(default: %(default)s)",
    )
    args = parser.parse_args(argv)

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not args.concordance.exists():
        print(
            f"[ERROR] concordance.csv not found: {args.concordance}\n"
            f"        Run build_concordance.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    low_vaf_labels = set(args.low_vaf_labels)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load and analyse
    # ------------------------------------------------------------------
    print(f"Loading concordance data from {args.concordance} ...")
    concordance_rows = load_concordance(args.concordance)
    print(f"  {len(concordance_rows)} rows loaded.")

    print("Extracting discordant variants ...")
    details = build_discordance_details(
        concordance_rows,
        depth_threshold=args.depth_threshold,
        low_vaf_labels=low_vaf_labels,
    )
    print(f"  {len(details)} discordant variants found.")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    details_out = args.out_dir / "discordance_details.csv"
    with details_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=DETAILS_FIELDS)
        writer.writeheader()
        writer.writerows(details)
    print(f"Discordance details written to {details_out}.")

    summary = build_discordance_summary(details)
    summary_out = args.out_dir / "discordance_summary.csv"
    with summary_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(summary)
    print(f"Discordance summary written to {summary_out}.")

    # ------------------------------------------------------------------
    # Stdout summary
    # ------------------------------------------------------------------
    print()
    print_summary(details, depth_threshold=args.depth_threshold)


if __name__ == "__main__":
    main()
