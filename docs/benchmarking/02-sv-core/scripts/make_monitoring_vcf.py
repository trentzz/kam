#!/usr/bin/env python3
"""Generate a monitoring-mode truth VCF from kam discovery results.

For each PASS SV call (LargeDeletion, TandemDuplication, Inversion) in the
discovery variants.tsv, this script applies the same variant-key extraction
logic as kam's targeting::extract_variant_key Rust function, producing a VCF
that monitoring mode can match exactly.

Usage:
    python make_monitoring_vcf.py <variants.tsv> <output.vcf>
"""

import sys


def complement(b: str) -> str:
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(b, b)


def is_reverse_complement(ref: str, alt: str) -> bool:
    if len(ref) < 2 or len(ref) != len(alt):
        return False
    return all(complement(r) == a for r, a in zip(ref, reversed(alt)))


def partial_inversion_len(ref: str, alt: str) -> int | None:
    """Detect the length of a central inverted segment, if one exists."""
    if len(ref) != len(alt):
        return None
    left = next((i for i, (r, a) in enumerate(zip(ref, alt)) if r != a), None)
    if left is None:
        return None
    right = len(ref) - 1 - next(
        i for i, (r, a) in enumerate(zip(reversed(ref), reversed(alt))) if r != a
    )
    if right < left:
        return None
    if is_reverse_complement(ref[left:right + 1], alt[left:right + 1]):
        return right - left + 1
    return None


def parse_target_id(target_id: str) -> tuple[str, int] | None:
    """Extract (chrom, 0-based target_start) from 'chrN:start-end[_label]'."""
    try:
        chrom, rest = target_id.split(':', 1)
        start_str = rest.split('-')[0]
        return chrom, int(start_str)
    except (ValueError, IndexError):
        return None


def extract_variant_key(
    target_id: str, ref_seq: str, alt_seq: str
) -> tuple[str, int, str, str] | None:
    """Reimplementation of kam's targeting::extract_variant_key.

    Returns (chrom, 1-based_pos, ref_allele, alt_allele) or None.
    """
    parsed = parse_target_id(target_id)
    if parsed is None:
        return None
    chrom, target_start = parsed

    if ref_seq == alt_seq:
        return None

    # Find leftmost differing position.
    diff_pos = next(
        (i for i, (r, a) in enumerate(zip(ref_seq, alt_seq)) if r != a), None
    )
    if diff_pos is None:
        return None

    ref_rem = ref_seq[diff_pos:]
    alt_rem = alt_seq[diff_pos:]
    max_cs = max(0, min(len(ref_rem), len(alt_rem)) - 1)
    common_suffix = 0
    for i in range(1, max_cs + 1):
        if ref_rem[-i] == alt_rem[-i]:
            common_suffix += 1
        else:
            break

    if common_suffix > 0:
        ref_trimmed = ref_seq[diff_pos: len(ref_seq) - common_suffix]
        alt_trimmed = alt_seq[diff_pos: len(alt_seq) - common_suffix]
    else:
        ref_trimmed = ref_seq[diff_pos:]
        alt_trimmed = alt_seq[diff_pos:]

    if len(ref_seq) == len(alt_seq):
        # SNV / MNV / Inversion — output the trimmed differing region.
        genomic_pos = target_start + diff_pos + 1  # 0-based offset → 1-based VCF
        return chrom, genomic_pos, ref_trimmed, alt_trimmed
    elif len(ref_seq) > len(alt_seq):
        return _deletion_key(chrom, target_start, diff_pos, ref_seq, ref_trimmed, alt_trimmed)
    else:
        return _insertion_key(chrom, target_start, diff_pos, ref_seq, alt_trimmed, ref_trimmed)


def _indel_key(chrom, target_start, diff_pos, ref_seq, alt_seq, ref_trimmed, alt_trimmed):
    """Route trimmed sequences to the appropriate indel handler."""
    # Remove inner common suffix.
    inner_cs = 0
    for i in range(1, min(len(ref_trimmed), len(alt_trimmed)) + 1):
        if ref_trimmed[-i] == alt_trimmed[-i]:
            inner_cs += 1
        else:
            break
    if inner_cs > 0:
        ref_min = ref_trimmed[:-inner_cs]
        alt_min = alt_trimmed[:-inner_cs]
    else:
        ref_min = ref_trimmed
        alt_min = alt_trimmed

    # Remove inner common prefix.
    inner_cp = 0
    for r, a in zip(ref_min, alt_min):
        if r == a:
            inner_cp += 1
        else:
            break
    ref_min = ref_min[inner_cp:]
    alt_min = alt_min[inner_cp:]
    indel_start = diff_pos + inner_cp

    if len(ref_seq) > len(alt_seq):
        return _deletion_key_from_min(chrom, target_start, indel_start, ref_seq, ref_min)
    else:
        return _insertion_key_from_min(chrom, target_start, indel_start, ref_seq, alt_min)


def _deletion_key(chrom, target_start, diff_pos, ref_seq, ref_trimmed, alt_trimmed):
    # Inner suffix removal.
    inner_cs = 0
    for i in range(1, min(len(ref_trimmed), len(alt_trimmed)) + 1):
        if ref_trimmed[-i] == alt_trimmed[-i]:
            inner_cs += 1
        else:
            break
    ref_min = ref_trimmed[:-inner_cs] if inner_cs else ref_trimmed
    alt_min = alt_trimmed[:-inner_cs] if inner_cs else alt_trimmed

    # Inner prefix removal.
    inner_cp = sum(1 for r, a in zip(ref_min, alt_min) if r == a)
    # Only prefix matches; stop as soon as mismatch found.
    inner_cp = 0
    for r, a in zip(ref_min, alt_min):
        if r == a:
            inner_cp += 1
        else:
            break
    ref_min = ref_min[inner_cp:]
    alt_min = alt_min[inner_cp:]
    indel_start = diff_pos + inner_cp
    return _deletion_key_from_min(chrom, target_start, indel_start, ref_seq, ref_min)


def _deletion_key_from_min(chrom, target_start, indel_start, ref_seq, del_seq_str):
    del_seq = list(del_seq_str)
    anchor_pos = indel_start - 1

    # Left-normalise: while last base of del_seq matches the base at anchor_pos.
    while anchor_pos > 0 and del_seq and ref_seq[anchor_pos] == del_seq[-1]:
        last = del_seq.pop()
        del_seq.insert(0, last)
        anchor_pos -= 1

    if anchor_pos >= 0:
        anchor = ref_seq[anchor_pos]
        ref_allele = anchor + ''.join(del_seq)
        alt_allele = anchor
        genomic_pos = target_start + anchor_pos + 1  # 1-based VCF
        return chrom, genomic_pos, ref_allele, alt_allele
    else:
        genomic_pos = target_start + indel_start + 1
        return chrom, genomic_pos, ''.join(del_seq), ''


def _insertion_key(chrom, target_start, diff_pos, ref_seq, ins_trimmed, ref_trimmed):
    # Inner suffix removal.
    inner_cs = 0
    for i in range(1, min(len(ref_trimmed), len(ins_trimmed)) + 1):
        if ref_trimmed[-i] == ins_trimmed[-i]:
            inner_cs += 1
        else:
            break
    ref_min = ref_trimmed[:-inner_cs] if inner_cs else ref_trimmed
    ins_min = ins_trimmed[:-inner_cs] if inner_cs else ins_trimmed

    inner_cp = 0
    for r, a in zip(ref_min, ins_min):
        if r == a:
            inner_cp += 1
        else:
            break
    ref_min = ref_min[inner_cp:]
    ins_min = ins_min[inner_cp:]
    indel_start = diff_pos + inner_cp
    return _insertion_key_from_min(chrom, target_start, indel_start, ref_seq, ins_min)


def _insertion_key_from_min(chrom, target_start, indel_start, ref_seq, ins_seq_str):
    ins_seq = list(ins_seq_str)
    anchor_pos = indel_start - 1

    while anchor_pos > 0 and ins_seq and ref_seq[anchor_pos] == ins_seq[-1]:
        last = ins_seq.pop()
        ins_seq.insert(0, last)
        anchor_pos -= 1

    if anchor_pos >= 0:
        anchor = ref_seq[anchor_pos]
        ref_allele = anchor
        alt_allele = anchor + ''.join(ins_seq)
        genomic_pos = target_start + anchor_pos + 1  # 1-based VCF
        return chrom, genomic_pos, ref_allele, alt_allele
    else:
        genomic_pos = target_start + indel_start + 1
        return chrom, genomic_pos, '', ''.join(ins_seq)


SV_TYPES = {'LargeDeletion', 'TandemDuplication', 'Inversion'}


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <variants.tsv> <output.vcf>", file=sys.stderr)
        sys.exit(1)

    tsv_path, vcf_path = sys.argv[1], sys.argv[2]
    entries = []

    with open(tsv_path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 14:
                continue
            target_id, variant_type, ref_seq, alt_seq = parts[0], parts[1], parts[2], parts[3]
            filter_col = parts[13]
            if filter_col != 'PASS':
                continue
            if variant_type not in SV_TYPES:
                continue
            key = extract_variant_key(target_id, ref_seq, alt_seq)
            if key is None:
                print(f"Warning: could not extract key for {target_id} {variant_type}", file=sys.stderr)
                continue
            chrom, pos, ref_allele, alt_allele = key
            entries.append((chrom, pos, ref_allele, alt_allele, variant_type))

    with open(vcf_path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n')
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for chrom, pos, ref_allele, alt_allele, sv_type in entries:
            svtype = {'LargeDeletion': 'DEL', 'TandemDuplication': 'DUP', 'Inversion': 'INV'}[sv_type]
            f.write(f'{chrom}\t{pos}\t.\t{ref_allele}\t{alt_allele}\t.\tPASS\tSVTYPE={svtype}\n')

    print(f"Wrote {len(entries)} SV monitoring entries to {vcf_path}")


if __name__ == '__main__':
    main()
