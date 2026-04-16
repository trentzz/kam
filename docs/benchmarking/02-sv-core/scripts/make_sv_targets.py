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
- Fusions/translocations (ALT contains BND bracket notation or SVTYPE=BND):
  junction target concatenates the partner A and partner B segments. The BND
  strand orientation (FF, FR, RF, RR) determines whether either segment is
  reverse-complemented before concatenation.

Requirements
------------
- pysam (for VCF parsing)
- pyfaidx or samtools faidx for reference access (optional, see --ref)

If pysam is not available, a minimal VCF parser is used.
"""

import argparse
import re
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
    # BND / fusion: recognised by SVTYPE=BND or bracket notation in ALT.
    if info.get("SVTYPE") == "BND" or "[" in alt or "]" in alt:
        return "BND"
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
    """Construct a reference-window target for a deletion.

    The target is the full genomic reference window spanning ``window // 2``
    bases of left context, the entire deleted region (with anchor), and
    ``window // 2`` bases of right context.  This is the reference sequence;
    the alt path (via junction k-mers from --sv-junctions) is shorter by
    exactly ``sv_len`` bases.

    Using the reference window as the target ensures that reference molecules
    are on-target (via the target k-mers in the allowlist) and that the
    reference path reconstructed by the DFS matches ``target_seq``, allowing
    the caller to correctly identify it as the reference allele.

    pos is 1-based (VCF convention).
    """
    half = window // 2
    # In 0-based coordinates:
    #   anchor   = pos - 1
    #   deleted  = [pos : pos + sv_len]   (100 bases for sv_len=100)
    #   right of deletion starts at pos + sv_len
    left_start = pos - 1 - half        # 50 bp before anchor
    right_end = pos + sv_len + half    # 50 bp after deletion end

    if left_start < 0:
        return None

    try:
        target_seq = fetch_reference(ref_path, chrom, left_start, right_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    if len(target_seq) < 2:
        return None

    target_id = f"{chrom}:{left_start}-{right_end}_DEL_{sv_len}bp"
    return target_id, target_seq


def make_duplication_junction(
    ref_path: Path, chrom: str, pos: int, sv_len: int, window: int
) -> tuple[str, str] | None:
    """Construct a reference-window target for a tandem duplication.

    The target is the genomic reference window spanning ``window // 2`` bases
    of left context, the duplicated region, and ``window // 2`` bases of right
    context.  The reference path matches this window.  The alt path (tandem
    duplicate) is longer by exactly ``sv_len`` bases and is found via the DUP
    junction k-mers in --sv-junctions.

    pos is 1-based (VCF convention).
    """
    half = window // 2
    dup_start = pos - 1       # 0-based start of duplicated region
    dup_end = pos - 1 + sv_len  # 0-based exclusive end
    left_start = dup_start - half
    right_end = dup_end + half

    if left_start < 0:
        return None

    try:
        target_seq = fetch_reference(ref_path, chrom, left_start, right_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    if len(target_seq) < 2:
        return None

    target_id = f"{chrom}:{left_start}-{right_end}_DUP_{sv_len}bp"
    return target_id, target_seq


def make_inversion_junction(
    ref_path: Path, chrom: str, pos: int, sv_len: int, window: int
) -> tuple[str, str] | None:
    """Construct a reference-window target for an inversion.

    The target is the genomic reference window spanning ``window // 2`` bases
    of left context, the inverted region, and ``window // 2`` bases of right
    context.  The reference path matches this window.  The alt path has the
    same length but the central region is replaced by the reverse complement,
    and it is found via the INV junction k-mers in --sv-junctions.

    pos is 1-based (VCF convention).
    """
    half = window // 2
    inv_start = pos - 1       # 0-based
    inv_end = pos - 1 + sv_len
    left_start = max(0, inv_start - half)
    right_end = inv_end + half

    try:
        target_seq = fetch_reference(ref_path, chrom, left_start, right_end)
    except Exception as e:
        print(f"  WARNING: could not fetch reference for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    if len(target_seq) < 2:
        return None

    target_id = f"{chrom}:{left_start}-{right_end}_INV_{sv_len}bp"
    return target_id, target_seq


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
    # Include anchor base in the left flank (left_end = pos, exclusive).
    left_end = pos
    left_start = left_end - flank
    # Right flank starts after deleted region: [pos:pos+sv_len] deleted → right at pos+sv_len.
    right_start = pos + sv_len
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


def make_inversion_junction_kmer_seq(
    ref_path: Path, chrom: str, pos: int, sv_len: int, k: int
) -> str | None:
    """Return synthetic junction sequences for an inversion (both breakpoints).

    An inversion at [inv_start, inv_end) creates junction k-mers at two sites:
      - Left breakpoint: last (k-1) bases of left flank + first (k-1) bases of
        rc(inverted region).
      - Right breakpoint: last (k-1) bases of rc(inverted region) + first (k-1)
        bases of right flank.

    Both junction sequences are concatenated and returned as a single FASTA
    entry.  They provide the alt-path edges the DFS needs to traverse the
    inversion without following the reference sequence.
    """
    flank = k - 1
    inv_start = pos - 1   # 0-based
    inv_end = pos - 1 + sv_len

    try:
        left_context = fetch_reference(ref_path, chrom, inv_start - flank, inv_start)
        inv_region = fetch_reference(ref_path, chrom, inv_start, inv_end)
        right_context = fetch_reference(ref_path, chrom, inv_end, inv_end + flank)
    except Exception as e:
        print(f"  WARNING: INV junction k-mer fetch failed for {chrom}:{pos}: {e}", file=sys.stderr)
        return None

    rc_inv = reverse_complement(inv_region)
    # Left breakpoint junction: last (k-1) of left context + first (k-1) of rc(inv).
    left_junc = left_context + rc_inv[:flank]
    # Right breakpoint junction: last (k-1) of rc(inv) + first (k-1) of right context.
    right_junc = rc_inv[-flank:] + right_context
    return left_junc + right_junc


# ── BND / fusion helpers ─────────────────────────────────────────────────────


# Regex to parse VCF BND ALT fields. The four forms are:
#   t]chr:pos]   (FF)
#   t[chr:pos[   (FR)
#   ]chr:pos]t   (RF)
#   [chr:pos[t   (RR)
_BND_RE = re.compile(
    r"^(?:"
    r"(?P<ff_base>[A-Za-z]+)\](?P<ff_chr>[^:]+):(?P<ff_pos>\d+)\]"
    r"|(?P<fr_base>[A-Za-z]+)\[(?P<fr_chr>[^:]+):(?P<fr_pos>\d+)\["
    r"|\](?P<rf_chr>[^:]+):(?P<rf_pos>\d+)\](?P<rf_base>[A-Za-z]+)"
    r"|\[(?P<rr_chr>[^:]+):(?P<rr_pos>\d+)\[(?P<rr_base>[A-Za-z]+)"
    r")$"
)


def parse_bnd_alt(alt: str) -> tuple[str, int, str] | None:
    """Parse a VCF BND ALT field and return (mate_chrom, mate_pos, orientation).

    mate_pos is 1-based (VCF convention). orientation is one of
    "FF", "FR", "RF", "RR".

    Returns None if the ALT field does not match any known BND pattern.
    """
    m = _BND_RE.match(alt)
    if m is None:
        return None
    if m.group("ff_chr"):
        return m.group("ff_chr"), int(m.group("ff_pos")), "FF"
    if m.group("fr_chr"):
        return m.group("fr_chr"), int(m.group("fr_pos")), "FR"
    if m.group("rf_chr"):
        return m.group("rf_chr"), int(m.group("rf_pos")), "RF"
    if m.group("rr_chr"):
        return m.group("rr_chr"), int(m.group("rr_pos")), "RR"
    return None


def make_fusion_junction(
    ref_path: Path,
    chrom_a: str,
    pos_a: int,
    chrom_b: str,
    pos_b: int,
    orientation: str,
    window: int,
) -> tuple[str, str] | None:
    """Construct a fusion junction target from two partner breakpoints.

    Each partner contributes ``window // 2`` bases around its breakpoint.
    Partner A contributes bases ending at the breakpoint; partner B contributes
    bases starting at the breakpoint.

    The orientation controls reverse-complementing of segments:
      - FF: both segments used as-is (forward-forward).
      - FR: partner B's segment is reverse-complemented.
      - RF: partner A's segment is reverse-complemented.
      - RR: both segments are reverse-complemented.

    pos_a and pos_b are 1-based (VCF convention).
    """
    half = window // 2

    # Partner A: half bases ending at the breakpoint.
    a_end = pos_a  # 0-based exclusive (pos_a in 0-based = pos_a - 1, plus 1 for exclusive)
    a_start = a_end - half
    if a_start < 0:
        return None

    # Partner B: half bases starting at the breakpoint.
    b_start = pos_b - 1  # 0-based inclusive
    b_end = b_start + half

    try:
        seg_a = fetch_reference(ref_path, chrom_a, a_start, a_end)
        seg_b = fetch_reference(ref_path, chrom_b, b_start, b_end)
    except Exception as e:
        print(
            f"  WARNING: could not fetch reference for BND "
            f"{chrom_a}:{pos_a} / {chrom_b}:{pos_b}: {e}",
            file=sys.stderr,
        )
        return None

    # Apply orientation-dependent reverse-complementing.
    if orientation == "FR":
        seg_b = reverse_complement(seg_b)
    elif orientation == "RF":
        seg_a = reverse_complement(seg_a)
    elif orientation == "RR":
        seg_a = reverse_complement(seg_a)
        seg_b = reverse_complement(seg_b)
    # FF: no transformation needed.

    junction_seq = seg_a + seg_b

    if len(junction_seq) < 2:
        return None

    orient_suffix = f"__{orientation}" if orientation != "FF" else ""
    target_id = (
        f"BND__{chrom_a}:{a_start}-{a_end}__{chrom_b}:{b_start}-{b_end}"
        f"__fusion{orient_suffix}"
    )
    return target_id, junction_seq


def make_fusion_junction_kmer_seq(
    ref_path: Path,
    chrom_a: str,
    pos_a: int,
    chrom_b: str,
    pos_b: int,
    orientation: str,
    k: int,
) -> str | None:
    """Return the synthetic junction k-mer sequence for a fusion breakpoint.

    This is ``k-1`` bases from partner A (ending at the breakpoint) joined with
    ``k-1`` bases from partner B (starting at the breakpoint). The orientation
    controls reverse-complementing exactly as in ``make_fusion_junction``.

    pos_a and pos_b are 1-based (VCF convention).
    """
    flank = k - 1

    a_end = pos_a
    a_start = a_end - flank
    if a_start < 0:
        return None

    b_start = pos_b - 1
    b_end = b_start + flank

    try:
        seg_a = fetch_reference(ref_path, chrom_a, a_start, a_end)
        seg_b = fetch_reference(ref_path, chrom_b, b_start, b_end)
    except Exception as e:
        print(
            f"  WARNING: BND junction k-mer fetch failed for "
            f"{chrom_a}:{pos_a} / {chrom_b}:{pos_b}: {e}",
            file=sys.stderr,
        )
        return None

    if orientation == "FR":
        seg_b = reverse_complement(seg_b)
    elif orientation == "RF":
        seg_a = reverse_complement(seg_a)
    elif orientation == "RR":
        seg_a = reverse_complement(seg_a)
        seg_b = reverse_complement(seg_b)

    return seg_a + seg_b


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

    # Track BND mate IDs already processed. Each BND pair produces two VCF
    # records; we only need to emit one junction target per pair.
    seen_bnd_mates: set[str] = set()

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

            # BND records have a separate code path: no SVLEN, partner
            # coordinates are parsed from the ALT field.
            if stype == "BND":
                parsed = parse_bnd_alt(alt)
                if parsed is None:
                    print(
                        f"  SKIP {chrom}:{pos} — could not parse BND ALT field: {alt}",
                        file=sys.stderr,
                    )
                    continue

                mate_chrom, mate_pos, orientation = parsed

                # Deduplicate BND mate pairs. Use a canonical key so that
                # whichever record we see first produces the junction.
                pair_key = tuple(sorted([
                    (chrom, pos),
                    (mate_chrom, mate_pos),
                ]))
                if pair_key in seen_bnd_mates:
                    continue
                seen_bnd_mates.add(pair_key)

                result = make_fusion_junction(
                    args.ref, chrom, pos, mate_chrom, mate_pos,
                    orientation, args.window,
                )
                if result is None:
                    print(
                        f"  SKIP {chrom}:{pos} — BND junction construction failed",
                        file=sys.stderr,
                    )
                    continue

                target_id, junction_window = result
                ft.write(f">{target_id}\n{junction_window}\n")
                n_written += 1

                junc_seq = make_fusion_junction_kmer_seq(
                    args.ref, chrom, pos, mate_chrom, mate_pos,
                    orientation, args.kmer_size,
                )
                if junc_seq is not None and len(junc_seq) >= args.kmer_size:
                    junc_id = f"{chrom}:{pos}_{mate_chrom}:{mate_pos}_BND_junction"
                    fj.write(f">{junc_id}\n{junc_seq}\n")
                    n_junction_kmers += 1

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
            elif stype == "INV":
                junc_seq = make_inversion_junction_kmer_seq(
                    args.ref, chrom, pos, length, args.kmer_size
                )

            if junc_seq is not None and len(junc_seq) >= args.kmer_size:
                junc_id = f"{chrom}:{pos}_{stype}_junction"
                fj.write(f">{junc_id}\n{junc_seq}\n")
                n_junction_kmers += 1

    print(f"Processed {n_total} SV records.")
    print(f"Wrote {n_written} junction targets to {args.out_targets}")
    print(f"Wrote {n_junction_kmers} junction k-mer sequences to {args.out_junctions}")


if __name__ == "__main__":
    main()
