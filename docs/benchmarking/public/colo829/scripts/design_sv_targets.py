#!/usr/bin/env python3
"""Design kam target and junction sequences for the COLO829 SV benchmark.

For each validated SV in the COLO829 truth VCF, extracts:
  1. A reference-flanking target sequence (100 bp either side of each
     breakpoint). Written to the targets FASTA for `kam run --targets`.
  2. A synthetic junction sequence that spans the two breakpoints with
     100 bp of each side fused. Written to the junctions FASTA for
     `kam run --sv-junctions`.

Junction sequences allow kam to index k-mers that span the structural
variant breakpoint. These k-mers are absent from the reference target
sequences but present in tumour reads that carry the rearrangement.

SV types handled:
  DEL  — single-chromosome deletion. Breakpoints are POS and END.
  DUP  — tandem duplication. Breakpoints are POS and END.
  INV  — inversion. Breakpoints are POS and END.
  INS  — insertion. Breakpoint is POS; inserted sequence from ALT field.
  TRA  — translocation. Two chromosomes; CHR2/END2 or MATEID in INFO.

Usage:
    python3 design_sv_targets.py \\
        --vcf  data/COLO829_truthset_somatic_v4.1.vcf \\
        --ref  /path/to/GRCh37.fa \\
        --targets-out  data/colo829_sv_targets.fa \\
        --junctions-out data/colo829_sv_junctions.fa

Arguments:
    --vcf           COLO829 truth VCF (plain, not gzipped by default).
    --ref           Reference genome FASTA (GRCh37, must have .fai index).
    --targets-out   Output FASTA path for reference target sequences.
    --junctions-out Output FASTA path for SV junction sequences.
    --flank         Flanking bases each side of each breakpoint (default: 100).

Notes:
    - Requires pysam: pip install pysam
    - Reference must be indexed: samtools faidx <ref.fa>
    - The COLO829 truth VCF uses GRCh37 chromosome names (1, 2, ...).
      Ensure your reference matches.
    - For translocations, only the left breakpoint junction is synthesised
      when CHR2/END2 information is absent from the INFO field.
"""

import argparse
import sys
from pathlib import Path


# SV types present in the COLO829 truth set.
SV_TYPES = {"DEL", "DUP", "INV", "INS", "TRA", "BND"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Design SV target and junction sequences for COLO829."
    )
    p.add_argument("--vcf", required=True, help="COLO829 truth VCF.")
    p.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA (GRCh37, indexed with samtools faidx).",
    )
    p.add_argument(
        "--targets-out",
        required=True,
        help="Output FASTA for reference-flanking target sequences.",
    )
    p.add_argument(
        "--junctions-out",
        required=True,
        help="Output FASTA for synthetic SV junction sequences.",
    )
    p.add_argument(
        "--flank",
        type=int,
        default=100,
        help="Flanking bases each side of each breakpoint (default: 100).",
    )
    return p.parse_args()


def parse_info(info_str: str) -> dict[str, str]:
    """Parse VCF INFO field into a dict."""
    result = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = "1"
    return result


def parse_vcf(vcf_path: str) -> list[dict]:
    """Parse the COLO829 truth VCF into a list of SV records."""
    records = []
    with open(vcf_path) as fh:
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
            info = parse_info(fields[7]) if len(fields) > 7 else {}

            sv_type = info.get("SVTYPE", "").upper()
            if not sv_type:
                # Infer from ALT field.
                if alt.startswith("<DEL"):
                    sv_type = "DEL"
                elif alt.startswith("<DUP"):
                    sv_type = "DUP"
                elif alt.startswith("<INV"):
                    sv_type = "INV"
                elif alt.startswith("<INS"):
                    sv_type = "INS"
                else:
                    sv_type = "UNK"

            # Parse END coordinate.
            end = pos
            if "END" in info:
                try:
                    end = int(info["END"])
                except ValueError:
                    pass
            elif "END2" in info:
                try:
                    end = int(info["END2"])
                except ValueError:
                    pass

            # For translocations, parse second chromosome from CHR2 or ALT.
            chrom2 = info.get("CHR2", chrom)
            if sv_type in ("TRA", "BND") and chrom2 == chrom:
                # Try to parse from BND ALT notation: T[chr2:pos[ or ]chr2:pos]T
                import re
                m = re.search(r"[[\]]([^[\]:]+):(\d+)[[\]]", alt)
                if m:
                    chrom2 = m.group(1)
                    end = int(m.group(2))

            records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "filter": filt,
                    "sv_type": sv_type,
                    "end": end,
                    "chrom2": chrom2,
                    "info": info,
                }
            )
    return records


def fetch_safe(ref, chrom: str, start: int, end: int) -> str | None:
    """Fetch reference sequence, trying chr-prefix normalisation.

    Returns uppercase sequence or None if the chromosome is not found.
    """
    if chrom not in ref.references:
        alt = chrom.lstrip("chr") if chrom.startswith("chr") else f"chr{chrom}"
        if alt in ref.references:
            chrom = alt
        else:
            return None
    length = ref.get_reference_length(chrom)
    start = max(0, start)
    end = min(end, length)
    if start >= end:
        return None
    return ref.fetch(chrom, start, end).upper()


def make_target_entries(
    records: list[dict], ref, flank: int
) -> list[tuple[str, str]]:
    """Extract reference-flanking target sequences for each SV breakpoint.

    For single-chromosome SVs: two windows around POS and END.
    For translocations: one window per chromosome.
    Returns (header, sequence) tuples.
    """
    entries = []
    seen: set[str] = set()

    def add_window(chrom: str, centre: int, label: str) -> None:
        start = centre - 1 - flank  # convert 1-based POS to 0-based start
        end = centre + flank
        key = f"{chrom}:{start}-{end}"
        if key in seen:
            return
        seen.add(key)
        seq = fetch_safe(ref, chrom, start, end)
        if seq:
            entries.append((f"{chrom}:{centre}__ref_flank__{label}", seq))

    for r in records:
        sv_type = r["sv_type"]
        add_window(r["chrom"], r["pos"], f"{sv_type}_left")
        if r["end"] != r["pos"] or r["chrom2"] != r["chrom"]:
            add_window(r["chrom2"], r["end"], f"{sv_type}_right")

    return entries


def make_junction_entries(
    records: list[dict], ref, flank: int
) -> list[tuple[str, str]]:
    """Synthesise junction sequences that span each SV breakpoint.

    For deletions and inversions: fuse flank bp of left breakpoint with
    flank bp of right breakpoint into a single sequence.
    For duplications: fuse end-flank with start+flank (tandem junction).
    For insertions: flanking sequence around insertion site.
    For translocations: fuse left breakpoint of chrom1 with right
    breakpoint of chrom2.

    Returns (header, sequence) tuples.
    """
    entries = []

    for i, r in enumerate(records):
        sv_type = r["sv_type"]
        chrom = r["chrom"]
        pos = r["pos"]
        end = r["end"]
        chrom2 = r["chrom2"]

        if sv_type in ("DEL", "INV"):
            # Left side: bases to the left of POS.
            left_seq = fetch_safe(ref, chrom, pos - 1 - flank, pos - 1)
            # Right side: bases to the right of END.
            right_seq = fetch_safe(ref, chrom, end, end + flank)
            if left_seq and right_seq:
                junction = left_seq + right_seq
                header = f"{chrom}:{pos}-{end}__{sv_type}_junction__{i}"
                entries.append((header, junction))

        elif sv_type == "DUP":
            # Tandem duplication junction: bases at END fused to bases at POS.
            left_seq = fetch_safe(ref, chrom, end - flank, end)
            right_seq = fetch_safe(ref, chrom, pos - 1, pos - 1 + flank)
            if left_seq and right_seq:
                junction = left_seq + right_seq
                header = f"{chrom}:{pos}-{end}__DUP_junction__{i}"
                entries.append((header, junction))

        elif sv_type == "INS":
            # For insertions with known sequence in ALT, embed it.
            alt = r["alt"]
            if alt.startswith("<") or len(alt) <= 1:
                # Symbolic ALT or unknown insertion sequence — skip junction.
                continue
            ins_seq = alt[1:]  # strip leading REF base
            left_seq = fetch_safe(ref, chrom, pos - 1 - flank, pos - 1)
            right_seq = fetch_safe(ref, chrom, pos, pos + flank)
            if left_seq and right_seq:
                junction = left_seq + ins_seq + right_seq
                header = f"{chrom}:{pos}__INS_junction__{i}"
                entries.append((header, junction))

        elif sv_type in ("TRA", "BND"):
            # Translocation: fuse left breakpoint of chrom with right of chrom2.
            left_seq = fetch_safe(ref, chrom, pos - 1 - flank, pos - 1)
            right_seq = fetch_safe(ref, chrom2, end, end + flank)
            if left_seq and right_seq:
                junction = left_seq + right_seq
                header = (
                    f"{chrom}:{pos}__{chrom2}:{end}__TRA_junction__{i}"
                )
                entries.append((header, junction))

    return entries


def write_fasta(entries: list[tuple[str, str]], path: str) -> int:
    """Write (header, sequence) pairs to a FASTA file."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for header, seq in entries:
            fh.write(f">{header}\n{seq}\n")
    return len(entries)


def main() -> None:
    args = parse_args()

    # Check inputs.
    if not Path(args.vcf).exists():
        print(f"[ERROR] VCF not found: {args.vcf}", file=sys.stderr)
        print(
            "        Download with: sh docs/benchmarking/public/scripts/download_colo829.sh",
            file=sys.stderr,
        )
        sys.exit(1)
    if not Path(args.ref).exists():
        print(f"[ERROR] Reference FASTA not found: {args.ref}", file=sys.stderr)
        sys.exit(1)
    fai = args.ref + ".fai"
    if not Path(fai).exists():
        print(f"[ERROR] Reference index not found: {fai}", file=sys.stderr)
        print(f"        Create it with: samtools faidx {args.ref}", file=sys.stderr)
        sys.exit(1)

    try:
        import pysam
    except ImportError:
        print("[ERROR] pysam is required. Install with: pip install pysam", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Parsing SV records from: {args.vcf}")
    records = parse_vcf(args.vcf)
    print(f"[INFO] Found {len(records)} SV records.")

    # Summarise by type.
    from collections import Counter
    type_counts = Counter(r["sv_type"] for r in records)
    for sv_type, count in sorted(type_counts.items()):
        print(f"  {sv_type}: {count}")

    ref = pysam.FastaFile(args.ref)

    print(f"[INFO] Designing reference-flanking targets (±{args.flank} bp)...")
    target_entries = make_target_entries(records, ref, args.flank)
    n_targets = write_fasta(target_entries, args.targets_out)
    print(f"[INFO] Written {n_targets} target sequences to: {args.targets_out}")

    print(f"[INFO] Synthesising SV junction sequences (±{args.flank} bp)...")
    junction_entries = make_junction_entries(records, ref, args.flank)
    n_junctions = write_fasta(junction_entries, args.junctions_out)
    print(f"[INFO] Written {n_junctions} junction sequences to: {args.junctions_out}")

    ref.close()

    if n_targets == 0:
        print(
            "[WARN] No target sequences written. Check chromosome naming "
            "conventions between VCF and reference.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("")
    print("Run kam with:")
    print(f"  kam run \\")
    print(f"      --targets {args.targets_out} \\")
    print(f"      --sv-junctions {args.junctions_out} \\")
    print(f"      ...")


if __name__ == "__main__":
    main()
