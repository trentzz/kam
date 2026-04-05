#!/usr/bin/env python3
"""Design kam target sequences for the SEQC2 HCC1395 benchmark.

For each PASS variant in the SEQC2 truth VCF, extracts a 100 bp window of
reference sequence centred on the variant position. The resulting FASTA is
suitable as input to `kam run --targets`.

Usage:
    python3 design_targets.py \\
        --vcf  data/HCC1395_truth_SNV.vcf.gz \\
        --ref  /path/to/GRCh38.fa \\
        --output data/seqc2_snv_targets.fa

    python3 design_targets.py \\
        --vcf  data/HCC1395_truth_indel.vcf.gz \\
        --ref  /path/to/GRCh38.fa \\
        --output data/seqc2_indel_targets.fa

Arguments:
    --vcf       SEQC2 truth VCF (gzipped or plain).
    --ref       Reference genome FASTA (must have .fai index; create with
                samtools faidx if absent).
    --output    Output FASTA path for kam target sequences.
    --flank     Flanking bases each side of variant (default: 100).
    --pass-only Consider only PASS variants in the truth VCF (default: true).

Notes:
    - The script requires pysam (pip install pysam).
    - The reference must be indexed: samtools faidx <ref.fa>
    - Overlapping windows are deduplicated by region string before fetching.
    - Chromosome name normalisation: the VCF and reference must use the same
      convention (chr1 vs 1). A mismatch warning is printed for each variant
      that cannot be fetched.
"""

import argparse
import gzip
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Design 100 bp flanking target sequences for SEQC2 truth variants."
    )
    p.add_argument("--vcf", required=True, help="SEQC2 truth VCF (gzip or plain).")
    p.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA (indexed with samtools faidx).",
    )
    p.add_argument("--output", required=True, help="Output FASTA path.")
    p.add_argument(
        "--flank",
        type=int,
        default=100,
        help="Flanking bases each side of the variant position (default: 100).",
    )
    p.add_argument(
        "--pass-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include only PASS variants (default: true).",
    )
    return p.parse_args()


def open_vcf(path: str):
    """Open a VCF file, handling gzip transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_vcf_variants(vcf_path: str, pass_only: bool) -> list[dict]:
    """Parse variants from a VCF file.

    Returns a list of dicts with keys: chrom, pos (1-based), ref, alt, filter.
    """
    variants = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ref = fields[3]
            alt = fields[4]
            filt = fields[6] if len(fields) > 6 else "."
            if pass_only and filt not in ("PASS", "."):
                continue
            variants.append(
                {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "filter": filt}
            )
    return variants


def design_targets(
    variants: list[dict],
    ref_path: str,
    flank: int,
) -> list[tuple[str, str]]:
    """Fetch flanking sequences from the reference for each variant.

    Returns a list of (header, sequence) tuples.

    Requires pysam. Each window extends [pos - flank, pos + flank] in 0-based
    half-open coordinates, clamped to chromosome bounds.
    """
    try:
        import pysam
    except ImportError:
        print(
            "[ERROR] pysam is required. Install with: pip install pysam",
            file=sys.stderr,
        )
        sys.exit(1)

    ref = pysam.FastaFile(ref_path)
    entries = []
    seen_regions: set[str] = set()
    skipped = 0

    for v in variants:
        chrom = v["chrom"]
        pos = v["pos"]  # 1-based
        # Convert to 0-based half-open for pysam.
        start = max(0, pos - 1 - flank)
        end = pos + flank  # pysam clamps to chrom length

        region_key = f"{chrom}:{start}-{end}"
        if region_key in seen_regions:
            continue
        seen_regions.add(region_key)

        # Verify chromosome exists.
        if chrom not in ref.references:
            # Try the opposite chr prefix convention.
            alt_chrom = chrom.lstrip("chr") if chrom.startswith("chr") else f"chr{chrom}"
            if alt_chrom in ref.references:
                chrom = alt_chrom
            else:
                skipped += 1
                if skipped <= 5:
                    print(
                        f"[WARN] Chromosome not in reference: {v['chrom']} "
                        f"(tried {alt_chrom}). Skipping.",
                        file=sys.stderr,
                    )
                continue

        seq = ref.fetch(chrom, start, end).upper()
        if not seq:
            skipped += 1
            continue

        header = (
            f"{chrom}:{start + 1}-{end}__"
            f"{v['ref']}>{v['alt']}__"
            f"pos{pos}"
        )
        entries.append((header, seq))

    if skipped:
        print(
            f"[WARN] Skipped {skipped} variant(s) due to missing chromosomes "
            f"or empty sequences.",
            file=sys.stderr,
        )

    ref.close()
    return entries


def write_fasta(entries: list[tuple[str, str]], output_path: str) -> int:
    """Write (header, sequence) pairs to a FASTA file.

    Returns the number of records written.
    """
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(f">{header}\n{seq}\n")
    return len(entries)


def main() -> None:
    args = parse_args()

    # Check inputs exist.
    if not Path(args.vcf).exists():
        print(f"[ERROR] VCF not found: {args.vcf}", file=sys.stderr)
        print(
            "        Download with: sh docs/benchmarking/05-public/scripts/download_seqc2.sh",
            file=sys.stderr,
        )
        sys.exit(1)
    if not Path(args.ref).exists():
        print(f"[ERROR] Reference FASTA not found: {args.ref}", file=sys.stderr)
        print(
            "        Provide the path to GRCh38.fa (or GRCh37.fa for GRCh37 coordinates).",
            file=sys.stderr,
        )
        sys.exit(1)
    fai = args.ref + ".fai"
    if not Path(fai).exists():
        print(f"[ERROR] Reference index not found: {fai}", file=sys.stderr)
        print(f"        Create it with: samtools faidx {args.ref}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Parsing variants from: {args.vcf}")
    variants = parse_vcf_variants(args.vcf, args.pass_only)
    print(f"[INFO] Found {len(variants)} {'PASS ' if args.pass_only else ''}variants.")

    if not variants:
        print("[ERROR] No variants found. Check --pass-only and VCF format.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Fetching {args.flank} bp flanking sequences from: {args.ref}")
    entries = design_targets(variants, args.ref, args.flank)

    n = write_fasta(entries, args.output)
    print(f"[INFO] Written {n} target sequences to: {args.output}")

    if n == 0:
        print(
            "[WARN] No targets written. The VCF and reference may use different "
            "chromosome naming conventions.",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
