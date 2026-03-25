"""
build_concordance.py

Join alignment-based results with kam results to produce a per-variant
concordance table.

Inputs
------
  alignment_baseline.csv
      One row per (sample, variant). Columns include alignment_detected,
      alignment_vaf, alignment_depth.

  kam_results/{sample_id}/variants.tsv
      Per-sample kam output in tumour-informed (monitoring) mode.
      Columns: target_id, variant_type, ref_seq, alt_seq, vaf, vaf_ci_low,
               vaf_ci_high, n_molecules_ref, n_molecules_alt, n_duplex_alt,
               n_simplex_alt, strand_bias_p, confidence, filter

Outputs
-------
  concordance.csv
      One row per (sample, variant). Columns:
        sample_id, input_ng, vaf_level, chromosome, position, ref, alt,
        variant_type, alignment_detected, alignment_vaf,
        kam_detected, kam_vaf, kam_confidence,
        concordance

      concordance values:
        both_detected    — alignment and kam both called the variant
        alignment_only   — only alignment detected it
        kam_only         — only kam detected it
        neither          — neither method detected it

  concordance_summary.csv
      Aggregate concordance rates broken down by vaf_level and variant_type.

Usage
-----
  python build_concordance.py [--baseline PATH] [--kam-dir PATH] [--out-dir PATH]

Defaults assume the script is run from anywhere; paths resolve relative to
the script's location inside docs/benchmarking/comparison/scripts/.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
COMPARISON_DIR = SCRIPT_DIR.parent

DEFAULT_BASELINE = COMPARISON_DIR / "alignment_baseline.csv"
DEFAULT_KAM_DIR = COMPARISON_DIR / "kam_results"
DEFAULT_OUT_DIR = COMPARISON_DIR


# ---------------------------------------------------------------------------
# VAF label ordering for summary tables
# ---------------------------------------------------------------------------

# Used to sort rows in summary output predictably.
VAF_ORDER = [
    "0pc",
    "0p001pc",
    "0p01pc",
    "0p1pc",
    "0p25pc",
    "0p5pc",
    "1pc",
    "2pc",
]


def _vaf_sort_key(label: str) -> float:
    """Convert VAF label to a numeric value for sorting."""
    vaf_map = {
        "0pc": 0.0,
        "0p001pc": 0.00001,
        "0p01pc": 0.0001,
        "0p1pc": 0.001,
        "0p25pc": 0.0025,
        "0p5pc": 0.005,
        "1pc": 0.01,
        "2pc": 0.02,
    }
    return vaf_map.get(label, float("inf"))


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

_TARGET_ID_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def _parse_target_id(target_id: str) -> tuple[str, int, int] | None:
    """Return (chromosome, start, end) from a target_id like chr1:100-200.

    Returns None if the target_id does not match the expected format.
    """
    m = _TARGET_ID_RE.match(target_id.strip())
    if not m:
        return None
    return m.group("chrom"), int(m.group("start")), int(m.group("end"))


def _position_in_target(pos: int, start: int, end: int) -> bool:
    """Return True if pos falls within [start, end] (1-based inclusive)."""
    return start <= pos <= end


def _parse_variant_position_from_seqs(
    target_id: str, ref_seq: str, alt_seq: str
) -> int | None:
    """Infer the variant position from the ref and alt sequences.

    The reference sequence for a target spans [start, end] (100 bp). The
    variant position is identified by finding where ref_seq and alt_seq first
    differ. The 1-based genomic position is then start + offset.

    Returns None if target_id cannot be parsed or the sequences do not differ.
    """
    parsed = _parse_target_id(target_id)
    if parsed is None:
        return None
    chrom, start, _end = parsed

    # Find first differing position.
    min_len = min(len(ref_seq), len(alt_seq))
    for i in range(min_len):
        if ref_seq[i] != alt_seq[i]:
            return start + i  # 1-based: start is already 1-based

    # Sequences share a common prefix equal to the shorter one — difference is
    # an insertion or deletion relative to ref.
    return start + min_len


def _extract_alleles_from_seqs(
    ref_seq: str, alt_seq: str, offset: int
) -> tuple[str, str]:
    """Extract REF and ALT alleles from full-length target sequences.

    offset is the 0-based index within the sequences where they first differ.
    For SNVs this returns single characters. For indels it returns the
    minimal allele representation (no common prefix beyond offset-1).
    """
    ref_part = ref_seq[offset:]
    alt_part = alt_seq[offset:]

    # Find common suffix to trim.
    i = 0
    while i < len(ref_part) and i < len(alt_part) and ref_part[-(i + 1)] == alt_part[-(i + 1)]:
        i += 1

    if i > 0:
        ref_allele = ref_part[: len(ref_part) - i]
        alt_allele = alt_part[: len(alt_part) - i]
    else:
        ref_allele = ref_part
        alt_allele = alt_part

    # Return at least one character (anchor base not needed here — the position
    # already pins us to the correct locus).
    return ref_allele or ref_seq[offset], alt_allele or alt_seq[offset]


# ---------------------------------------------------------------------------
# Load alignment baseline
# ---------------------------------------------------------------------------

def load_alignment_baseline(path: Path) -> dict[tuple[str, str, int, str, str], dict]:
    """Load alignment_baseline.csv into a lookup dict.

    Key: (sample_id, chromosome, position, ref, alt)
    Value: full row dict.
    """
    data: dict[tuple[str, str, int, str, str], dict] = {}

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            key = (
                row["sample_id"],
                row["chromosome"],
                int(row["position"]),
                row["ref"],
                row["alt"],
            )
            data[key] = row

    return data


# ---------------------------------------------------------------------------
# Load kam results for one sample
# ---------------------------------------------------------------------------

def load_kam_results(tsv_path: Path) -> list[dict]:
    """Load a kam variants.tsv file.

    Returns a list of dicts with keys:
      chromosome, position, ref, alt, vaf, confidence, filter
    """
    results = []
    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            target_id = row.get("target_id", "").strip()
            ref_seq = row.get("ref_seq", "").strip()
            alt_seq = row.get("alt_seq", "").strip()

            if not target_id or not ref_seq or not alt_seq:
                continue

            parsed = _parse_target_id(target_id)
            if parsed is None:
                continue
            chrom, start, _end = parsed

            # Find offset where sequences differ.
            min_len = min(len(ref_seq), len(alt_seq))
            offset = None
            for i in range(min_len):
                if ref_seq[i] != alt_seq[i]:
                    offset = i
                    break
            if offset is None:
                offset = min_len

            pos = start + offset  # 1-based genomic position

            ref_allele, alt_allele = _extract_alleles_from_seqs(ref_seq, alt_seq, offset)

            try:
                vaf = float(row.get("vaf", ""))
            except (ValueError, TypeError):
                vaf = None

            results.append({
                "chromosome": chrom,
                "position": pos,
                "ref": ref_allele,
                "alt": alt_allele,
                "vaf": vaf,
                "confidence": row.get("confidence", ""),
                "filter": row.get("filter", ""),
            })

    return results


# ---------------------------------------------------------------------------
# Build concordance table
# ---------------------------------------------------------------------------

def build_concordance(
    baseline: dict[tuple[str, str, int, str, str], dict],
    kam_dir: Path,
    sample_ids: list[str],
) -> list[dict]:
    """Join alignment baseline with kam results for each sample.

    For each (sample, variant) in the baseline, find the matching kam call (if
    any) and classify concordance.
    """
    rows = []

    for sample_id in sample_ids:
        # Load kam results for this sample.
        tsv_path = kam_dir / sample_id / "variants.tsv"
        if not tsv_path.exists():
            # kam has not been run for this sample. All rows get kam_detected=False.
            kam_calls: dict[tuple[str, int, str, str], dict] = {}
        else:
            raw = load_kam_results(tsv_path)
            # Index by (chromosome, position, ref, alt).
            kam_calls = {
                (r["chromosome"], r["position"], r["ref"], r["alt"]): r
                for r in raw
            }

        # Find all baseline rows for this sample.
        sample_rows = [
            row for key, row in baseline.items()
            if key[0] == sample_id
        ]

        for base_row in sample_rows:
            chrom = base_row["chromosome"]
            pos = int(base_row["position"])
            ref = base_row["ref"]
            alt = base_row["alt"]

            # Alignment result.
            align_detected_str = base_row.get("alignment_detected", "False")
            align_detected = align_detected_str.strip().upper() in ("TRUE", "1", "YES")
            align_vaf_raw = base_row.get("alignment_vaf", "")
            try:
                align_vaf = float(align_vaf_raw) if align_vaf_raw.strip() else None
            except ValueError:
                align_vaf = None

            # kam result.
            kam_key = (chrom, pos, ref, alt)
            kam_row = kam_calls.get(kam_key)
            if kam_row is not None:
                kam_detected = True
                kam_vaf = kam_row["vaf"]
                kam_confidence = kam_row["confidence"]
            else:
                kam_detected = False
                kam_vaf = None
                kam_confidence = ""

            # Concordance category.
            if align_detected and kam_detected:
                concordance = "both_detected"
            elif align_detected and not kam_detected:
                concordance = "alignment_only"
            elif not align_detected and kam_detected:
                concordance = "kam_only"
            else:
                concordance = "neither"

            rows.append({
                "sample_id": sample_id,
                "input_ng": base_row["input_ng"],
                "vaf_level": base_row["vaf_level"],
                "chromosome": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": base_row["variant_type"],
                "alignment_detected": align_detected,
                "alignment_vaf": align_vaf if align_vaf is not None else "",
                "kam_detected": kam_detected,
                "kam_vaf": kam_vaf if kam_vaf is not None else "",
                "kam_confidence": kam_confidence,
                "concordance": concordance,
            })

    # Deterministic sort: sample_id, then chromosome, then position.
    rows.sort(key=lambda r: (r["sample_id"], r["chromosome"], int(r["position"])))

    return rows


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def build_summary(rows: list[dict]) -> list[dict]:
    """Compute per-VAF and per-variant-type concordance rates.

    Returns a list of dicts suitable for writing to concordance_summary.csv.
    """
    # Aggregate by (vaf_level, variant_type).
    counts: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {"total": 0, "both_detected": 0, "alignment_only": 0, "kam_only": 0, "neither": 0}
    )

    for r in rows:
        key = (r["vaf_level"], r["variant_type"])
        counts[key]["total"] += 1
        counts[key][r["concordance"]] += 1

    summary_rows = []
    for (vaf_level, vtype), c in sorted(
        counts.items(),
        key=lambda x: (_vaf_sort_key(x[0][0]), x[0][1]),
    ):
        total = c["total"]
        if total == 0:
            continue

        align_detected = c["both_detected"] + c["alignment_only"]
        kam_detected = c["both_detected"] + c["kam_only"]

        summary_rows.append({
            "vaf_level": vaf_level,
            "variant_type": vtype,
            "total_calls": total,
            "both_detected": c["both_detected"],
            "alignment_only": c["alignment_only"],
            "kam_only": c["kam_only"],
            "neither": c["neither"],
            "alignment_sensitivity": (
                f"{align_detected / total:.4f}" if total > 0 else ""
            ),
            "kam_sensitivity": (
                f"{kam_detected / total:.4f}" if total > 0 else ""
            ),
            "concordance_rate": (
                f"{c['both_detected'] / max(align_detected, 1):.4f}"
                if align_detected > 0
                else ""
            ),
        })

    return summary_rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

CONCORDANCE_FIELDS = [
    "sample_id",
    "input_ng",
    "vaf_level",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant_type",
    "alignment_detected",
    "alignment_vaf",
    "kam_detected",
    "kam_vaf",
    "kam_confidence",
    "concordance",
]

SUMMARY_FIELDS = [
    "vaf_level",
    "variant_type",
    "total_calls",
    "both_detected",
    "alignment_only",
    "kam_only",
    "neither",
    "alignment_sensitivity",
    "kam_sensitivity",
    "concordance_rate",
]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build per-variant concordance table from alignment and kam results."
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=DEFAULT_BASELINE,
        help="Path to alignment_baseline.csv (default: %(default)s)",
    )
    parser.add_argument(
        "--kam-dir",
        type=Path,
        default=DEFAULT_KAM_DIR,
        help="Directory containing per-sample kam results (default: %(default)s)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for concordance CSV files (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    if not args.baseline.exists():
        print(
            f"[ERROR] Alignment baseline not found: {args.baseline}",
            file=sys.stderr,
        )
        sys.exit(1)

    if not args.kam_dir.exists():
        print(
            f"[ERROR] kam results directory not found: {args.kam_dir}\n"
            f"        Run run_kam_titration.sh first to generate kam results.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load baseline
    # ------------------------------------------------------------------
    print(f"Loading alignment baseline from {args.baseline} ...")
    baseline = load_alignment_baseline(args.baseline)
    print(f"  {len(baseline)} (sample, variant) rows loaded.")

    # Determine set of unique sample IDs from the baseline.
    sample_ids = sorted({key[0] for key in baseline})
    print(f"  {len(sample_ids)} unique samples.")

    # Warn about samples with no kam output.
    missing_samples = [
        sid for sid in sample_ids
        if not (args.kam_dir / sid / "variants.tsv").exists()
    ]
    if missing_samples:
        print(
            f"\n[WARN] kam results not found for {len(missing_samples)} sample(s):",
            file=sys.stderr,
        )
        for sid in missing_samples:
            print(f"       {sid}", file=sys.stderr)
        print(
            "       Run run_kam_titration.sh to generate missing results.\n"
            "       These samples will show kam_detected=False for all variants.",
            file=sys.stderr,
        )

    # ------------------------------------------------------------------
    # Build concordance table
    # ------------------------------------------------------------------
    print("\nBuilding concordance table ...")
    concordance_rows = build_concordance(baseline, args.kam_dir, sample_ids)
    print(f"  {len(concordance_rows)} rows built.")

    concordance_out = args.out_dir / "concordance.csv"
    with concordance_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CONCORDANCE_FIELDS)
        writer.writeheader()
        writer.writerows(concordance_rows)
    print(f"Concordance table written to {concordance_out}.")

    # ------------------------------------------------------------------
    # Build summary statistics
    # ------------------------------------------------------------------
    print("\nBuilding summary statistics ...")
    summary_rows = build_summary(concordance_rows)

    summary_out = args.out_dir / "concordance_summary.csv"
    with summary_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Summary written to {summary_out}.")

    # Print a quick overall summary to stdout.
    total = len(concordance_rows)
    if total > 0:
        counts = defaultdict(int)
        for r in concordance_rows:
            counts[r["concordance"]] += 1

        align_total = counts["both_detected"] + counts["alignment_only"]
        kam_total = counts["both_detected"] + counts["kam_only"]

        print(f"\nOverall ({total} calls):")
        print(f"  both_detected   : {counts['both_detected']:>5}  ({100 * counts['both_detected'] / total:.1f}%)")
        print(f"  alignment_only  : {counts['alignment_only']:>5}  ({100 * counts['alignment_only'] / total:.1f}%)")
        print(f"  kam_only        : {counts['kam_only']:>5}  ({100 * counts['kam_only'] / total:.1f}%)")
        print(f"  neither         : {counts['neither']:>5}  ({100 * counts['neither'] / total:.1f}%)")
        print(f"  alignment sensitivity: {align_total}/{total} = {100 * align_total / total:.1f}%")
        print(f"  kam sensitivity      : {kam_total}/{total} = {100 * kam_total / total:.1f}%")


if __name__ == "__main__":
    main()
