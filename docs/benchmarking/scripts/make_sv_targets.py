#!/usr/bin/env python3
"""Generate SV junction target FASTA and junction k-mer allowlist FASTA.

For each structural variant in a truth VCF, this script:
  1. Extracts a 100 bp window centred on the breakpoint junction from the
     reference genome.
  2. Writes the junction window to a target FASTA for use with ``--targets``.
  3. Optionally writes a synthetic junction sequence for each SV to a second
     FASTA for use with ``--sv-junctions``.

The junction k-mers for SVs (e.g. deletion alt junctions, duplication junctions)
are not in the reference sequence and are therefore not in the standard allowlist.
The ``--sv-junctions`` FASTA lets ``kam index`` augment the allowlist with these
k-mers so that SV alt paths are not invisible to the de Bruijn graph.

Usage
-----
    python make_sv_targets.py \\
        --vcf truth_svs.vcf \\
        --ref GRCh38.fa \\
        --out-targets sv_targets.fa \\
        --out-junctions sv_junctions.fa \\
        [--window 100] \\
        [--kmer-size 31]

Input VCF
---------
The truth VCF must have standard fields (CHROM, POS, REF, ALT) and may use
either full-sequence alleles or symbolic alleles (<DEL>, <DUP>, <INV>).
For symbolic alleles, SVLEN in the INFO field is required to determine the
breakpoint coordinates.

Supported SV types
------------------
- Deletions (ALT = <DEL> or alt shorter than ref): junction target spans the
  deletion breakpoint (left flank + right flank).
- Tandem duplications (ALT = <DUP> or alt longer than ref): junction target
  spans the duplication boundary (end of duplicated region + start).
- Inversions (ALT = <INV> or same length): junction target spans the
  inversion breakpoint (left flank + rev-comp of the inversion start).

Requirements
------------
- pysam (for VCF parsing)
- pyfaidx or samtools faidx for reference access (optional, see --ref)

If pysam is not available, a minimal VCF parser is used.
"""

import argparse
import sys
from pathlib import Path


# ── Constants ─────────────────────────────────────────────────────────────────

DEFAULT_WINDOW = 100
DEFAULT_KMER_SIZE = 31
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


# ── Helpers ────────────────────────────────────────────────────────────────────


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def fetch_reference(ref_path: Path, chrom: str, start: int, end: int) -> str:
    """Fetch reference sequence using pyfaidx or samtools faidx.

    start is 0-based, end is exclusive.
    """
    try:
        import pyfaidx  # type: ignore

        fa = pyfaidx.Fasta(str(ref_path))
        seq = fa[chrom][start:end].seq.upper()
        return seq
    except ImportError:
        pass

    # Fall back to samtools faidx.
    import subprocess

    region = f"{chrom}:{start + 1}-{end}"  # samtools is 1-based
    result = subprocess.run(
        ["samtools", "faidx", str(ref_path), region],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools faidx failed for {region}: {result.stderr.strip()}"
        )
    lines = result.stdout.strip().split("\n")
    return "".join(lines[1:]).upper()


def parse_vcf(vcf_path: Path):
    """Parse a VCF file and yield (chrom, pos, ref, alt, info_dict) tuples.

    pos is 1-based (VCF convention).
    """
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_str, _, ref, alt = parts[:5]
            info_str = parts[7] if len(parts) > 7 else "."
            pos = int(pos_str)
            info: dict[str, str] = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = "1"
            yield chrom, pos, ref, alt, info


def sv_type(ref: str, alt: str, info: dict) -> str | None:
    """Determine SV type from REF, ALT, and INFO fields."""
    if alt.upper() in ("<DEL>", "DEL") or info.get("SVTYPE") == "DEL":
        return "DEL"
    if alt.upper() in ("<DUP>", "DUP") or info.get("SVTYPE") == "DUP":
        return "DUP"
    if alt.upper() in ("<INV>", "INV") or info.get("SVTYPE") == "INV":
        return "INV"
    # Sequence-level classification.
    if len(ref) > len(alt) + 49:  # deletion ≥ 50 bp
        return "DEL"
    if len(alt) > len(ref) + 49:  # duplication ≥ 50 bp
        return "DUP"
    if len(ref) == len(alt) and len(ref) >= 50:
        return "INV"
    return None


def svlen(ref: str, alt: str, info: dict, stype: str) -> int:
    """Return the absolute SV length in bp."""
    if "SVLEN" in info:
        return abs(int(info["SVLEN"]))
    if stype == "DEL":
        return len(ref) - len(alt)
    if stype == "DUP":
        return len(alt) - len(ref)
    if stype == "INV":
        return len(ref)
    return 0


# ── Junction target construction ──────────────────────────────────────────────


def make_deletion_junction(
    ref_path: Path, chrom: str, pos: int, sv_len: int, window: int
) -> tuple[str, str] | None:
    """Construct a junction window for a deletion.

    The junction is the point where the left flank meets the right flank after
    the deleted sequence is removed.  The window captures ``window // 2`` bp of
    left flank and ``window // 2`` bp of right flank.

    pos is 1-based (VCF convention); the deleted region starts at pos.
    """
    half = window // 2
    # Left flank: bases immediately before the deleted region (0-based: pos-1-half to pos-1).
    left_start = pos - 1 - half
    left_end = pos - 1  # exclusive
    # Right flank: bases immediately after the deleted region.
    right_start = pos - 1 + sv_len
    right_end = right_start + half

    if left_start < 0:
        return None

    try:
        left = fetch_reference(ref_path, chrom, left_start, left_end)
        right = fetch_reference(ref_path, chrom, right_start, right_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    junction_seq = left + right
    if len(junction_seq) < 2:
        return None

    target_id = f"{chrom}:{left_start}-{right_end}_DEL_{sv_len}bp_maxpath{window + sv_len}"
    return target_id, junction_seq


def make_duplication_junction(
    ref_path: Path, chrom: str, pos: int, sv_len: int, window: int
) -> tuple[str, str] | None:
    """Construct a junction window for a tandem duplication.

    The junction k-mers are at the boundary where the duplicated segment meets
    a copy of itself.  The target captures the last ``half`` bp of the
    duplicated region + the first ``half`` bp of the duplicated region again.
    """
    half = window // 2
    # Duplicated region: pos-1 to pos-1+sv_len (0-based, exclusive end).
    dup_start = pos - 1
    dup_end = pos - 1 + sv_len

    # Junction: last half of dup + first half of dup.
    seg_end_start = max(dup_start, dup_end - half)
    seg_end_end = dup_end
    seg_begin_start = dup_start
    seg_begin_end = dup_start + half

    try:
        seg_end = fetch_reference(ref_path, chrom, seg_end_start, seg_end_end)
        seg_begin = fetch_reference(ref_path, chrom, seg_begin_start, seg_begin_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    junction_seq = seg_end + seg_begin
    if len(junction_seq) < 2:
        return None

    target_id = f"{chrom}:{dup_start}-{dup_end}_DUP_{sv_len}bp"
    return target_id, junction_seq


def make_inversion_junction(
    ref_path: Path, chrom: str, pos: int, sv_len: int, window: int
) -> tuple[str, str] | None:
    """Construct a junction window for an inversion.

    The left breakpoint junction captures the left flank + the rev-comp of the
    first ``half`` bp of the inverted region.
    """
    half = window // 2
    inv_start = pos - 1  # 0-based
    inv_end = pos - 1 + sv_len

    left_flank_start = max(0, inv_start - half)
    left_flank_end = inv_start
    inv_region_start = inv_start
    inv_region_end = min(inv_end, inv_start + half)

    try:
        left_flank = fetch_reference(ref_path, chrom, left_flank_start, left_flank_end)
        inv_region = fetch_reference(ref_path, chrom, inv_region_start, inv_region_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    junction_seq = left_flank + reverse_complement(inv_region)
    if len(junction_seq) < 2:
        return None

    target_id = f"{chrom}:{left_flank_start}-{inv_end}_INV_{sv_len}bp"
    return target_id, junction_seq


# ── Junction k-mer allowlist FASTA ────────────────────────────────────────────


def make_deletion_junction_kmer_seq(
    ref_path: Path, chrom: str, pos: int, sv_len: int, k: int
) -> str | None:
    """Return the synthetic junction sequence for a deletion.

    This is the concatenation of ``k-1`` bases before the deletion with
    ``k-1`` bases after the deletion.  The k-mers spanning this junction are
    the alt-specific k-mers for the deletion.
    """
    flank = k - 1
    left_end = pos - 1
    left_start = left_end - flank
    right_start = pos - 1 + sv_len
    right_end = right_start + flank

    if left_start < 0:
        return None

    try:
        left = fetch_reference(ref_path, chrom, left_start, left_end)
        right = fetch_reference(ref_path, chrom, right_start, right_end)
    except Exception as e:
        print(f"  WARNING: junction k-mer fetch failed for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    return left + right


def make_duplication_junction_kmer_seq(
    ref_path: Path, chrom: str, pos: int, sv_len: int, k: int
) -> str | None:
    """Return the synthetic junction sequence for a tandem duplication.

    The duplication junction k-mers span the end of the duplicated segment
    joined to its own beginning: last (k-1) bp of dup + first (k-1) bp of dup.
    """
    flank = k - 1
    dup_start = pos - 1
    dup_end = pos - 1 + sv_len

    seg_end_start = max(dup_start, dup_end - flank)
    seg_end_end = dup_end
    seg_begin_start = dup_start
    seg_begin_end = min(dup_end, dup_start + flank)

    try:
        seg_end = fetch_reference(ref_path, chrom, seg_end_start, seg_end_end)
        seg_begin = fetch_reference(ref_path, chrom, seg_begin_start, seg_begin_end)
    except Exception as e:
        print(f"  WARNING: junction k-mer fetch failed for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    return seg_end + seg_begin


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate SV junction targets and junction k-mer allowlist FASTA files"
    )
    parser.add_argument("--vcf", required=True, type=Path, help="Truth SV VCF file")
    parser.add_argument(
        "--ref", required=True, type=Path, help="Reference genome FASTA (indexed)"
    )
    parser.add_argument(
        "--out-targets",
        required=True,
        type=Path,
        help="Output target FASTA for --targets",
    )
    parser.add_argument(
        "--out-junctions",
        required=True,
        type=Path,
        help="Output junction k-mer FASTA for --sv-junctions",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help=f"Junction window size in bp (default: {DEFAULT_WINDOW})",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=DEFAULT_KMER_SIZE,
        help=f"K-mer size (default: {DEFAULT_KMER_SIZE})",
    )
    args = parser.parse_args()

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    if not args.ref.exists():
        print(f"ERROR: reference not found: {args.ref}", file=sys.stderr)
        sys.exit(1)

    n_total = 0
    n_written = 0
    n_junction_kmers = 0

    with open(args.out_targets, "w") as ft, open(args.out_junctions, "w") as fj:
        for chrom, pos, ref, alt, info in parse_vcf(args.vcf):
            n_total += 1
            stype = sv_type(ref, alt, info)
            if stype is None:
                print(
                    f"  SKIP {chrom}:{pos} — not a supported SV type (ref={ref[:10]}, alt={alt[:10]})",
                    file=sys.stderr,
                )
                continue

            length = svlen(ref, alt, info, stype)
            if length < 1:
                print(f"  SKIP {chrom}:{pos} — could not determine SVLEN", file=sys.stderr)
                continue

            # Build junction target window.
            result: tuple[str, str] | None = None
            if stype == "DEL":
                result = make_deletion_junction(args.ref, chrom, pos, length, args.window)
            elif stype == "DUP":
                result = make_duplication_junction(args.ref, chrom, pos, length, args.window)
            elif stype == "INV":
                result = make_inversion_junction(args.ref, chrom, pos, length, args.window)

            if result is None:
                print(f"  SKIP {chrom}:{pos} — junction window construction failed", file=sys.stderr)
                continue

            target_id, junction_window = result
            ft.write(f">{target_id}\n{junction_window}\n")
            n_written += 1

            # Build junction k-mer sequence for the allowlist.
            junc_seq: str | None = None
            if stype == "DEL":
                junc_seq = make_deletion_junction_kmer_seq(
                    args.ref, chrom, pos, length, args.kmer_size
                )
            elif stype == "DUP":
                junc_seq = make_duplication_junction_kmer_seq(
                    args.ref, chrom, pos, length, args.kmer_size
                )
            # For inversions, the junction window itself serves as the allowlist entry.
            elif stype == "INV":
                junc_seq = junction_window

            if junc_seq is not None and len(junc_seq) >= args.kmer_size:
                junc_id = f"{chrom}:{pos}_{stype}_junction"
                fj.write(f">{junc_id}\n{junc_seq}\n")
                n_junction_kmers += 1

    print(f"Processed {n_total} SV records.")
    print(f"Wrote {n_written} junction targets to {args.out_targets}")
    print(f"Wrote {n_junction_kmers} junction k-mer sequences to {args.out_junctions}")


if __name__ == "__main__":
    main()
