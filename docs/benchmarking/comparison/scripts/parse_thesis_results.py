"""
parse_thesis_results.py

Parse thesis alignment-based variant detection results into a standardised
format for comparison with kam output.

Data provenance
---------------
Source: /mnt/tzeng-local/tzeng-thesis/

Key inputs:
  - benchmark-dataframes/dedup-v2-generate-v3/   (SNV per-sample TSV files)
  - benchmark-dataframes-indel/dedup-v2/          (indel per-sample TSV files)
  - titration.probes.QC.pass.tsv                  (truth variant set, 375 rows)

Pipeline chosen: dedup-v2-generate-v3 (SNV) and dedup-v2 (indel).
These correspond to the HUMID v2 deduplication + RASCALL v3 calling pipeline,
which is the primary alignment-based result used in the thesis.

File format (per-sample TSV)
----------------------------
Columns present when a variant is detected:
  Location                  chr:pos (e.g. chr1:43349338)
  Tumour fraction           VAF label (e.g. 2pc, 0p1pc)
  Input amount (ng)         float (e.g. 15, 30)
  Variant type (benchmark)  SNP or INDEL
  Found variant             TRUE / FALSE
  Target sequence length    70bp
  RASCALL rVAF              float, detected VAF (NA when not detected)
  RASCALL Min_coverage      float, depth at locus (NA when not detected)
  ... (other RASCALL columns, not used)

When a variant is not detected, columns beyond "Target sequence length" are
empty (tab-separated but blank).

Outputs
-------
  alignment_baseline.csv   one row per (sample, variant)
  truth_variants.csv       one row per truth variant

Usage
-----
  python parse_thesis_results.py [--data-root PATH] [--out-dir PATH]

Defaults assume the script is run from anywhere; paths resolve relative to
the fixed thesis data directory.
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

THESIS_ROOT = Path("/mnt/tzeng-local/tzeng-thesis")

SNV_DIR = THESIS_ROOT / "benchmark-dataframes" / "dedup-v2-generate-v3"
INDEL_DIR = THESIS_ROOT / "benchmark-dataframes-indel" / "dedup-v2"
TRUTH_FILE = THESIS_ROOT / "titration.probes.QC.pass.tsv"

DEFAULT_OUT_DIR = Path(__file__).resolve().parent.parent  # docs/benchmarking/comparison/

# Map filename VAF tokens to numeric float values (as actual fraction, not percent).
VAF_MAP = {
    "0pc": 0.0,
    "0p001pc": 0.00001,
    "0p01pc": 0.0001,
    "0p1pc": 0.001,
    "0p25pc": 0.0025,
    "0p5pc": 0.005,
    "1pc": 0.01,
    "2pc": 0.02,
}


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

# Example: TWIST_STDV2_15ng_VAF_2pc_DEDUPED_70bp-targets_results.tsv
_FILENAME_RE = re.compile(
    r"TWIST_STDV2_(?P<ng>\d+)ng_VAF_(?P<vaf>[^_]+)_"
)


def parse_filename(name: str) -> tuple[int, float, str]:
    """Return (input_ng, vaf_fraction, vaf_label) from a results filename.

    Raises ValueError if the filename does not match the expected pattern.
    """
    m = _FILENAME_RE.search(name)
    if not m:
        raise ValueError(f"Cannot parse filename: {name!r}")
    ng = int(m.group("ng"))
    vaf_label = m.group("vaf")
    vaf = VAF_MAP.get(vaf_label)
    if vaf is None:
        raise ValueError(f"Unknown VAF token {vaf_label!r} in {name!r}")
    return ng, vaf, vaf_label


# ---------------------------------------------------------------------------
# Truth variant loading
# ---------------------------------------------------------------------------

def load_truth_variants(path: Path) -> dict[tuple[str, int], dict]:
    """Load truth variants keyed by (chromosome, position).

    The file uses quoted fields and has columns:
      Number, Chromosome, Position, Ref, Alt, Type

    Returns a dict mapping (chrom, pos) -> {chromosome, position, ref, alt, variant_type}.
    """
    truth = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quotechar='"')
        for row in reader:
            chrom = row["Chromosome"].strip('"')
            pos = int(str(row["Position"]).strip('"'))
            ref = row["Ref"].strip('"')
            alt = row["Alt"].strip('"')
            vtype = row["Type"].strip('"')
            # Normalise type label: the per-sample files use SNP/INDEL.
            # The truth file uses SNP/INDEL (same convention).
            truth[(chrom, pos)] = {
                "chromosome": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": vtype,
            }
    return truth


# ---------------------------------------------------------------------------
# Per-sample TSV parsing
# ---------------------------------------------------------------------------

def _parse_bool(value: str) -> bool:
    """Convert TRUE/FALSE/True/False string to bool."""
    return value.strip().upper() == "TRUE"


def _parse_float_or_na(value: str):
    """Return float, or None if the value is empty or not numeric."""
    v = value.strip()
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parse_sample_file(
    path: Path,
    input_ng: int,
    vaf_fraction: float,
    vaf_label: str,
    truth: dict[tuple[str, int], dict],
) -> list[dict]:
    """Parse one per-sample results TSV and return a list of row dicts.

    Each row corresponds to one truth variant in this sample.
    Variants not present in the truth set are skipped (should not occur in
    these benchmark files, but we guard against it).

    Columns in output:
      sample_id, input_ng, vaf_level, chromosome, position, ref, alt,
      variant_type, alignment_detected, alignment_vaf, alignment_depth
    """
    sample_id = path.stem  # filename without extension
    rows = []

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            location = row.get("Location", "").strip()
            if not location:
                continue

            # Location is chr:pos
            parts = location.split(":")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            # Look up Ref/Alt from truth set.
            key = (chrom, pos)
            truth_entry = truth.get(key)
            if truth_entry is None:
                # Position not in truth — skip silently.
                continue

            detected = _parse_bool(row.get("Found variant", "FALSE"))

            # RASCALL rVAF is the alignment-based VAF estimate.
            # RASCALL Min_coverage is the depth (minimum coverage across the
            # target sequence), used as a proxy for alignment depth.
            alignment_vaf = _parse_float_or_na(row.get("RASCALL rVAF", ""))
            alignment_depth = _parse_float_or_na(row.get("RASCALL Min_coverage", ""))

            rows.append({
                "sample_id": sample_id,
                "input_ng": input_ng,
                "vaf_level": vaf_label,
                "chromosome": chrom,
                "position": pos,
                "ref": truth_entry["ref"],
                "alt": truth_entry["alt"],
                "variant_type": truth_entry["variant_type"],
                "alignment_detected": detected,
                "alignment_vaf": alignment_vaf if alignment_vaf is not None else "",
                "alignment_depth": (
                    int(round(alignment_depth))
                    if alignment_depth is not None
                    else ""
                ),
            })

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

OUTPUT_FIELDNAMES = [
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
    "alignment_depth",
]

TRUTH_FIELDNAMES = [
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant_type",
]


def collect_results(
    data_dir: Path,
    truth: dict[tuple[str, int], dict],
) -> list[dict]:
    """Collect all per-sample results from a benchmark-dataframes directory."""
    all_rows: list[dict] = []
    tsv_files = sorted(data_dir.glob("*.tsv"))
    if not tsv_files:
        print(f"  WARNING: no TSV files found in {data_dir}", file=sys.stderr)
        return all_rows

    for tsv in tsv_files:
        try:
            ng, vaf_frac, vaf_label = parse_filename(tsv.name)
        except ValueError as exc:
            print(f"  SKIP {tsv.name}: {exc}", file=sys.stderr)
            continue

        rows = parse_sample_file(tsv, ng, vaf_frac, vaf_label, truth)
        all_rows.extend(rows)
        print(f"  {tsv.name}: {len(rows)} variants")

    return all_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    parser.add_argument(
        "--data-root",
        type=Path,
        default=THESIS_ROOT,
        help="Root directory of thesis data (default: %(default)s)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for CSV files (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    thesis_root: Path = args.data_root
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    snv_dir = thesis_root / "benchmark-dataframes" / "dedup-v2-generate-v3"
    indel_dir = thesis_root / "benchmark-dataframes-indel" / "dedup-v2"
    truth_file = thesis_root / "titration.probes.QC.pass.tsv"

    # ------------------------------------------------------------------
    # 1. Load truth variants
    # ------------------------------------------------------------------
    print(f"Loading truth variants from {truth_file} ...")
    truth = load_truth_variants(truth_file)
    print(f"  {len(truth)} truth variants loaded.")

    # Write truth_variants.csv
    truth_out = out_dir / "truth_variants.csv"
    with truth_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TRUTH_FIELDNAMES)
        writer.writeheader()
        for entry in sorted(truth.values(), key=lambda r: (r["chromosome"], r["position"])):
            writer.writerow({k: entry[k] for k in TRUTH_FIELDNAMES})
    print(f"Truth variants written to {truth_out} ({len(truth)} rows).")

    # ------------------------------------------------------------------
    # 2. Parse SNV results
    # ------------------------------------------------------------------
    print(f"\nParsing SNV results from {snv_dir} ...")
    snv_rows = collect_results(snv_dir, truth)
    print(f"  Total SNV rows: {len(snv_rows)}")

    # ------------------------------------------------------------------
    # 3. Parse indel results
    # ------------------------------------------------------------------
    print(f"\nParsing indel results from {indel_dir} ...")
    indel_rows = collect_results(indel_dir, truth)
    print(f"  Total indel rows: {len(indel_rows)}")

    # ------------------------------------------------------------------
    # 4. Merge and write alignment_baseline.csv
    # ------------------------------------------------------------------
    all_rows = snv_rows + indel_rows

    # Sort for deterministic output: chromosome order is lexicographic here.
    # Position is numeric. Then by input_ng, vaf_level.
    all_rows.sort(key=lambda r: (
        r["chromosome"],
        r["position"],
        r["input_ng"],
        r["vaf_level"],
    ))

    baseline_out = out_dir / "alignment_baseline.csv"
    with baseline_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nAlignment baseline written to {baseline_out} ({len(all_rows)} rows).")

    # Summary
    detected = sum(1 for r in all_rows if r["alignment_detected"])
    total = len(all_rows)
    print(f"Detection summary: {detected}/{total} ({100 * detected / total:.1f}%) detected.")


if __name__ == "__main__":
    main()
