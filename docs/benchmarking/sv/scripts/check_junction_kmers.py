#!/usr/bin/env python3
"""Diagnostic: check how many reads contain each junction k-mer from the targets.

This helps understand if the alt junction k-mers are present in the FASTQs.
"""

import gzip
import sys
from pathlib import Path

K = 31


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def canonical(seq: str) -> str:
    rc = rev_comp(seq)
    return seq if seq < rc else rc


def read_fasta(path: str) -> dict:
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


def get_junction_kmers(seq: str, k: int = K) -> set:
    """Return the k-mers that span the left/right boundary (alt-specific)."""
    # Junction is left_flank + right_flank. Half = len/2.
    half = len(seq) // 2
    # Junction k-mers start at half - k + 1 and end at half.
    start = max(0, half - k + 1)
    end = half
    kmers = set()
    for i in range(start, end + 1):
        if i + k <= len(seq):
            km = seq[i:i+k]
            if "N" not in km.upper():
                kmers.add(canonical(km.upper()))
    return kmers


def count_reads_with_kmer(fastq_gz: str, kmer_set: set, k: int = K) -> dict:
    counts = {km: 0 for km in kmer_set}
    with gzip.open(fastq_gz, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip()
            f.readline()  # +
            f.readline()  # qual
            for i in range(len(seq) - k + 1):
                km = seq[i:i+k].upper()
                if "N" not in km:
                    ckm = canonical(km)
                    if ckm in counts:
                        counts[ckm] += 1
    return counts


def main():
    targets_fa = "docs/benchmarking/sv/data/sv_targets.fa"
    r1_050 = "docs/benchmarking/sv/results/sim_vaf050/SV_VAF050_R1.fastq.gz"

    targets = read_fasta(targets_fa)
    for tid, seq in targets.items():
        jkmers = get_junction_kmers(seq)
        print(f"\nTarget: {tid} ({len(seq)} bp)")
        print(f"  Junction k-mers ({len(jkmers)} unique):")
        counts = count_reads_with_kmer(r1_050, jkmers)
        for km, cnt in sorted(counts.items(), key=lambda x: -x[1])[:5]:
            print(f"    {km[:20]}... : {cnt} reads")
        total = sum(counts.values())
        present = sum(1 for v in counts.values() if v > 0)
        print(f"  Total junction k-mer hits: {total}, k-mers present: {present}/{len(jkmers)}")


if __name__ == "__main__":
    main()
