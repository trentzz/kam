#!/usr/bin/env python3
"""Generate a synthetic reference FASTA and truth VCF for SV benchmarking.

Creates a 2000 bp reference on chr1 with known sequence composition designed
to support unambiguous SV detection:
  - Positions 100–299: region where a 100 bp deletion is introduced
  - Positions 500–699: region where a 100 bp tandem duplication is introduced
  - Positions 900–999: region where a 100 bp inversion is introduced

The sequence is designed with enough complexity (not ACGT repeat) that k-mers
spanning the breakpoints are unique.
"""

import random
import sys

SEED = 42
REF_LEN = 2000
CHROM = "chr1"


def make_complex_seq(length: int, rng: random.Random) -> str:
    """Generate a random sequence with no trivial repeats."""
    bases = "ACGT"
    seq = []
    for _ in range(length):
        seq.append(rng.choice(bases))
    return "".join(seq)


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def write_fasta(path: str, chrom: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{chrom}\n")
        # Write in 60-char lines.
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")


def write_fai(fasta_path: str, chrom: str, seq: str) -> None:
    """Write a samtools-compatible .fai index."""
    # Offset: length of ">chrom\n"
    offset = len(f">{chrom}\n")
    line_bases = 60
    line_width = 61  # 60 bases + newline
    length = len(seq)
    with open(fasta_path + ".fai", "w") as f:
        f.write(f"{chrom}\t{length}\t{offset}\t{line_bases}\t{line_width}\n")


def write_vcf(path: str, chrom: str, svs: list) -> None:
    """Write a truth VCF with full-sequence alleles for each SV.

    svs: list of (pos_1based, ref_seq, alt_seq, sv_type, sv_len) tuples.
    """
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={chrom},length={REF_LEN}>\n")
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        f.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref_seq, alt_seq, sv_type, sv_len, vaf in svs:
            svlen_str = -sv_len if sv_type == "DEL" else sv_len
            info = f"SVTYPE={sv_type};SVLEN={svlen_str};VAF={vaf}"
            f.write(f"{chrom}\t{pos}\t.\t{ref_seq}\t{alt_seq}\t.\tPASS\t{info}\n")


def main() -> None:
    rng = random.Random(SEED)
    ref = make_complex_seq(REF_LEN, rng)

    # ── Deletion: 100 bp deleted starting at pos 200 (1-based) ──────────────
    # VCF REF = anchor base + deleted sequence, ALT = anchor base only.
    del_pos = 200        # 1-based
    del_len = 100
    del_anchor = ref[del_pos - 1]             # base at pos 200
    del_ref = ref[del_pos - 1 : del_pos - 1 + del_len + 1]  # anchor + 100 deleted bases
    del_alt = del_anchor

    # ── Tandem duplication: 100 bp duplicated starting at pos 500 ────────────
    # VCF REF = anchor base, ALT = anchor base + duplicated sequence.
    dup_pos = 500        # 1-based
    dup_len = 100
    dup_anchor = ref[dup_pos - 1]
    dup_region = ref[dup_pos - 1 : dup_pos - 1 + dup_len]  # the duplicated sequence
    dup_ref = dup_anchor
    dup_alt = dup_anchor + dup_region  # insert a copy

    # ── Inversion: 100 bp inverted starting at pos 900 ────────────────────────
    # Full-sequence: REF = the 100 bp, ALT = reverse complement.
    inv_pos = 900        # 1-based
    inv_len = 100
    inv_ref = ref[inv_pos - 1 : inv_pos - 1 + inv_len]
    inv_alt = reverse_complement(inv_ref)

    out_dir = "docs/benchmarking/sv/data"

    ref_path = f"{out_dir}/ref.fa"
    write_fasta(ref_path, CHROM, ref)
    write_fai(ref_path, CHROM, ref)
    print(f"Wrote reference: {ref_path} ({REF_LEN} bp)")

    # Write truth VCFs at different VAFs.
    for vaf in [0.005, 0.01, 0.02, 0.05]:
        svs = [
            (del_pos, del_ref, del_alt, "DEL", del_len, vaf),
            (dup_pos, dup_ref, dup_alt, "DUP", dup_len, vaf),
            (inv_pos, inv_ref, inv_alt, "INV", inv_len, vaf),
        ]
        vcf_path = f"{out_dir}/truth_svs_vaf{int(vaf*1000):03d}.vcf"
        write_vcf(vcf_path, CHROM, svs)
        print(f"Wrote truth VCF: {vcf_path}")


if __name__ == "__main__":
    main()
