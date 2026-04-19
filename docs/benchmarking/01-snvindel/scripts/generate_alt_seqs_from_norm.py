#!/usr/bin/env python3
"""Generate alt-allele sequences from normalised targets and normalised VCF.

Replicates build_alt_seq logic from kam/src/rescue.rs so that the alt k-mers
added to the allowlist via --alt-as-ref exactly match what the rescue probe
queries at runtime.

For each variant in the VCF, finds the overlapping target window in the target
FASTA and applies the variant to produce an alt sequence. Writes the alt
sequences to a FASTA file.

Target headers are expected in the form "chrN:start-end" where start uses the
same coordinate convention as multiseqex (target_start = VCF_POS - flank,
which is effectively 1-based-like, matching the formula used by build_alt_seq:
offset = vcf_pos_1based - target_start).
"""

import sys
import re
import argparse
from pathlib import Path


def read_fasta(path):
    """Return list of (header, seq) pairs."""
    entries = []
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            entries.append((header, "".join(seq_parts)))
    return entries


def parse_target_header(header):
    """Parse 'chrN:start-end [extra]' → (chrom, start_int, end_int).

    start uses the multiseqex convention (vcf_pos - flank), matching the
    build_alt_seq formula: offset = vcf_pos_1based - start.
    """
    m = re.match(r"^(chr\S+):(\d+)-(\d+)", header)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def read_vcf_variants(path):
    """Return list of (chrom, pos_1based, ref, alt) tuples."""
    variants = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            variants.append((chrom, pos, ref, alt))
    return variants


def build_alt_seq(ref_seq, target_start, vcf_pos_1based, vcf_ref, vcf_alt):
    """Replicate build_alt_seq from kam/src/rescue.rs.

    target_start uses the multiseqex convention (vcf_pos - flank), so:
        offset = vcf_pos_1based - target_start

    Returns alt sequence bytes or None if ref check fails or variant outside window.
    """
    offset = vcf_pos_1based - target_start
    if offset < 0:
        return None
    if offset + len(vcf_ref) > len(ref_seq):
        return None
    ref_window = ref_seq[offset:offset + len(vcf_ref)]
    if ref_window.upper() != vcf_ref.upper():
        return None
    return ref_seq[:offset] + vcf_alt + ref_seq[offset + len(vcf_ref):]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", required=True, type=Path,
                        help="Target FASTA (e.g. targets_norm_100bp.fa)")
    parser.add_argument("--vcf", required=True, type=Path,
                        help="Normalised truth VCF")
    parser.add_argument("--output", required=True, type=Path,
                        help="Output alt-sequence FASTA")
    args = parser.parse_args()

    targets = read_fasta(args.targets)
    variants = read_vcf_variants(args.vcf)

    # Build a map from (chrom, start) → (header, seq) for fast lookup.
    target_map = {}
    for header, seq in targets:
        parsed = parse_target_header(header)
        if parsed is None:
            print(f"WARNING: cannot parse target header: {header}", file=sys.stderr)
            continue
        chrom, start, end = parsed
        target_map.setdefault(chrom, []).append((start, end, header, seq))

    n_written = 0
    n_failed = 0
    seen = set()  # deduplicate (chrom, start, vcf_ref, vcf_alt)

    with open(args.output, "w") as out:
        for chrom, pos, vcf_ref, vcf_alt in variants:
            candidates = target_map.get(chrom, [])
            matched = False
            for start, end, header, ref_seq in candidates:
                # Check if variant falls within this target window.
                # offset = pos - start; variant at ref_seq[offset..offset+len(vcf_ref)]
                offset = pos - start
                if offset < 0 or offset + len(vcf_ref) > len(ref_seq):
                    continue
                key = (chrom, start, vcf_ref, vcf_alt)
                if key in seen:
                    matched = True
                    break
                alt_seq = build_alt_seq(ref_seq, start, pos, vcf_ref, vcf_alt)
                if alt_seq is None:
                    n_failed += 1
                    print(
                        f"  FAIL ref check: {chrom}:{pos} {vcf_ref}>{vcf_alt} "
                        f"in target {header} (offset={offset}, "
                        f"window={ref_seq[offset:offset+len(vcf_ref)]})",
                        file=sys.stderr,
                    )
                    continue
                seen.add(key)
                out_header = f">{chrom}:{start}-{end} REF={vcf_ref} ALT={vcf_alt}"
                out.write(out_header + "\n")
                # Write 60 bases per line.
                for i in range(0, len(alt_seq), 60):
                    out.write(alt_seq[i:i+60] + "\n")
                n_written += 1
                matched = True
                break
            if not matched:
                print(
                    f"  NO TARGET: {chrom}:{pos} {vcf_ref}>{vcf_alt}",
                    file=sys.stderr,
                )
                n_failed += 1

    print(
        f"Written {n_written} alt sequences, {n_failed} failures.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
