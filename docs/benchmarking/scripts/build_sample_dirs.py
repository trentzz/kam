#!/usr/bin/env python3
"""Build per-sample result directories for benchmark suites.

For each benchmarked sample, creates a directory under
docs/benchmarking/{suite}/samples/{type}_vaf{tag}_{rep}/ containing:

  config.yaml           — varforge config used to generate the sample
  varforge_cmd.txt      — exact varforge command to regenerate the sample
  truth.tsv             — truth mutations (tab-separated)
  discovery.tsv         — kam discovery calls (tab-separated)
  tumour_informed.tsv   — kam tumour-informed calls (tab-separated)

Usage:
  python3 docs/benchmarking/build_sample_dirs.py

Run from the repository root. Skips samples where the source files do not
exist (e.g. where kam has not been run yet).
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent

VAF_TAGS = [
    "0005", "0010", "0015", "0020", "0025", "0030", "0035", "0040",
    "0050", "0060", "0075", "0100", "0125", "0150", "0175", "0200",
    "0250", "0300", "0350", "0400", "0500", "0600", "0700", "0800", "1000",
]
REPS = ["a", "b"]

# ─── VCF helpers ──────────────────────────────────────────────────────────────

def parse_info(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO field into a dict."""
    result: dict[str, str] = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = ""
    return result


def vcf_to_truth_tsv(vcf_path: Path, out_path: Path) -> None:
    """Convert a truth VCF to TSV.

    Output columns: chrom, pos, ref, alt, vaf, type
    """
    rows: list[str] = ["chrom\tpos\tref\talt\tvaf\ttype"]
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parse_info(parts[7]) if len(parts) > 7 else {}
            vaf = info.get("VAF", "")
            variant_type = info.get("SVTYPE", info.get("TYPE", ""))
            rows.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vaf}\t{variant_type}")
    out_path.write_text("\n".join(rows) + "\n")


def vcf_to_calls_tsv(vcf_path: Path, out_path: Path) -> None:
    """Convert a kam calls VCF to TSV.

    Output columns:
      chrom, pos, ref, alt, filter,
      vaf, vaf_lo, vaf_hi, nref, nalt, ndupalt, nsimalt, sbp, conf
    """
    header = (
        "chrom\tpos\tref\talt\tfilter\t"
        "vaf\tvaf_lo\tvaf_hi\tnref\tnalt\tndupalt\tnsimalt\tsbp\tconf"
    )
    rows: list[str] = [header]
    with vcf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos = parts[0], parts[1]
            ref, alt = parts[3], parts[4]
            filt = parts[6]
            info = parse_info(parts[7])
            vaf     = info.get("VAF", "")
            vaf_lo  = info.get("VAF_LO", "")
            vaf_hi  = info.get("VAF_HI", "")
            nref    = info.get("NREF", "")
            nalt    = info.get("NALT", "")
            ndupalt = info.get("NDUPALT", "")
            nsimalt = info.get("NSIMALT", "")
            sbp     = info.get("SBP", "")
            conf    = info.get("CONF", "")
            rows.append(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{filt}\t"
                f"{vaf}\t{vaf_lo}\t{vaf_hi}\t{nref}\t{nalt}\t"
                f"{ndupalt}\t{nsimalt}\t{sbp}\t{conf}"
            )
    out_path.write_text("\n".join(rows) + "\n")


# ─── Suite definitions ────────────────────────────────────────────────────────

def sv_samples() -> list[tuple[str, str, str, str]]:
    """Return (suite, type, tag, rep) tuples for the SV suite."""
    out = []
    for t in ["sv", "ins", "invdel"]:
        for tag in VAF_TAGS:
            for rep in REPS:
                out.append(("sv", t, tag, rep))
    return out


def snvindel_samples() -> list[tuple[str, str, str, str]]:
    """Return (suite, type, tag, rep) tuples for the SNV/indel suite."""
    out = []
    for t in ["snv", "indel"]:
        for tag in VAF_TAGS:
            for rep in REPS:
                out.append(("snvindel", t, tag, rep))
    return out


def truth_vcf_path(suite: str, vtype: str, tag: str, rep: str) -> Path:
    """Locate the truth VCF in the data directory."""
    data = REPO / "docs" / "benchmarking" / suite / "data"
    # SV suite uses type-specific prefixes.
    if suite == "sv":
        prefixes = {
            "sv":     "truth_svs",
            "ins":    "truth_ins",
            "invdel": "truth_invdel",
        }
        prefix = prefixes[vtype]
    else:
        prefixes = {
            "snv":   "truth_snvs",
            "indel": "truth_indels",
        }
        prefix = prefixes[vtype]
    return data / f"{prefix}_vaf{tag}_{rep}.vcf"


def build_sample(suite: str, vtype: str, tag: str, rep: str) -> bool:
    """Build one per-sample directory. Returns True on success."""
    base = REPO / "docs" / "benchmarking" / suite

    config_src = base / "configs" / f"{vtype}_vaf{tag}_{rep}.yaml"
    truth_src  = truth_vcf_path(suite, vtype, tag, rep)
    disc_src   = base / "results" / f"kam_{vtype}_vaf{tag}_{rep}" / "calls_discovery.vcf"
    ti_src     = base / "results" / f"kam_{vtype}_vaf{tag}_{rep}" / "calls_tumour_informed.vcf"

    # Skip if source files are absent (sample not yet generated or run).
    missing = [p for p in [config_src, truth_src, disc_src, ti_src] if not p.exists()]
    if missing:
        return False

    sample_dir = base / "samples" / f"{vtype}_vaf{tag}_{rep}"
    sample_dir.mkdir(parents=True, exist_ok=True)

    # config.yaml
    shutil.copy2(config_src, sample_dir / "config.yaml")

    # varforge_cmd.txt
    config_rel = config_src.relative_to(REPO)
    (sample_dir / "varforge_cmd.txt").write_text(
        f"# Run from the repository root.\n"
        f"varforge run {config_rel}\n"
    )

    # truth.tsv
    vcf_to_truth_tsv(truth_src, sample_dir / "truth.tsv")

    # discovery.tsv
    vcf_to_calls_tsv(disc_src, sample_dir / "discovery.tsv")

    # tumour_informed.tsv
    vcf_to_calls_tsv(ti_src, sample_dir / "tumour_informed.tsv")

    return True


def main() -> None:
    all_samples = sv_samples() + snvindel_samples()
    built = 0
    skipped = 0
    for suite, vtype, tag, rep in all_samples:
        if build_sample(suite, vtype, tag, rep):
            built += 1
        else:
            skipped += 1

    print(f"Built: {built}  Skipped (missing source files): {skipped}")


if __name__ == "__main__":
    main()
