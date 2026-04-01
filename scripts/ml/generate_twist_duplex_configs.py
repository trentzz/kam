#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for the Twist duplex ML dataset.

Produces 11,000 unique samples (10,000 train + 1,000 test) with:
- Fixed Twist UMI chemistry (length=5, duplex=true, inline=true)
- Randomised coverage, family size, PCR cycles, fragment size, and VAF
- One unique truth VCF per sample (1-3 variants each)

All outputs go to docs/benchmarking/ml-twist-duplex/.

Usage:
    python3 scripts/ml/generate_twist_duplex_configs.py [--dry-run]

Run from the repository root.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent

# ─── Output paths ─────────────────────────────────────────────────────────────

ML_DIR = REPO / "docs" / "benchmarking" / "ml-twist-duplex"
CONFIGS_TRAIN = ML_DIR / "configs" / "train"
CONFIGS_TEST = ML_DIR / "configs" / "test"
TRUTH_VCFS_DIR = ML_DIR / "truth_vcfs"

# ─── Reference files ──────────────────────────────────────────────────────────

SNVINDEL_REF = "docs/benchmarking/snvindel/data/ref.fa"
SV_REF = "docs/benchmarking/sv/data/ref.fa"

# Both references are chr1 with length 2000 (confirmed from ref.fa.fai).
# Variants are placed within safe margins from the ends.
CHROM = "chr1"
CHROM_LEN = 2000
# Minimum distance from contig edge to place a variant.
MARGIN = 200
# Maximum position for variant placement (leaving room for long alleles).
MAX_POS = CHROM_LEN - MARGIN

# Full reference sequence for snvindel ref (used when generating real REF bases).
SNVINDEL_REF_SEQ = (
    "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA"
    "CCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCG"
    "CAGTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTT"
    "TAAATGACCCTCTCGTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTAC"
    "AATGGCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAA"
    "ACAGCTCAGGAGCCAGTCCCCTACGTCGCATATCCTGGCCACTGGAGGTGAAGCGAATGG"
    "TATCGATACGTAGGAGGTGTGCCTTCGTAGGCTGTTTCTCAGGACGCCCAACTATTCTTT"
    "CCAATCCTACATCTGTTTCTTGCGTCGTAGCGGGACCCTCCATTGTTACTTATTAGGTTC"
    "TCGTTATGTCTCATAATCTCAGTGCTGGTGTGATAAGCAAACCACCCTACTGGCACGAAG"
    "TTCACAGAAGTGAGATTATGTCTCGTTTGGCAGTCTTGATGCTCGGGGGACACTTCTTTA"
    "AGCTCGGTGTGGTGGGCACGACCCTGGACGCGCGACGAAGCTAAGTTTGCAGTAATTAAC"
    "CGACATCTTTGTGAACCGACCCACATTTGACGGTACGCTACCGCAACGGTATGTGTTAAT"
    "GGAACAGACTTGCTTATGTGGACGTTGTATAGGGATATTACGTTACGCGTTAACCGATAC"
    "ATACTGGTTTCTCTCCAGTGGAGGTCTTGGTTGCCTCTAGTTTCTACGATATACTCATGG"
    "TAGTGTAACGCATAATCGAAGAGGGTCCTCCCATCTCCTGTGATGCATGGTGTGCTTACT"
    "GGGATGAATGCGCCGCAAGTAGCAGGTCCCGGCGTGGATACCTGATAGATGGTGACTAGC"
    "ATGTACAAGTAACCTTGTCTATTGAGCTTCGAGGATGCATACAAGCCCACCCGCAGCCGC"
    "AACAGCGACGACTAATTGATCAGTAATTTATTAAGCACGGTGTTAACTTCTGTTTAGTGG"
    "GCTAAAATAGCAGATGTAGGGACCTCAGGAGCTAGACGGGGACCTACAACTTTGCGGGAA"
    "CCAAGTTTTTGCAGTAGTGACTAACGCCGGGAATTCCTCGATATATAGTTTGATAGCTGA"
    "TACTTATGGCGCAACGGCCACGCCCACTTTGGCTATTGGAGAGTTAAGGAATTATCGTCA"
    "TAGACACTTCGGGTTGAGAGATGGCGACGGTCAGTGCATGAGGCCGTCCCCAGAAGCTCC"
    "CCTATGCTGTCCGTCGTTGTTCCCGATGAAGACGTCTACTGATATGCTAGCAGAGCCAGT"
    "CTTAAAGCCTAGCGAACTTAATACCGTAGCTCAGAATTATGGAGAGCAGCAGGCTTCCAT"
    "AGCACAGGTTGACGGAGGAGTTTTGCTTGGATATCGGAAGGGTTCTGTAGTGAATGCACT"
    "ACACGGTACTGGTACGTGGCAACTTAGGTCGTCACATCTAGGAGGCCGCACCCTAGGTCA"
    "AGTTTTACGATTGCCCTAACGCCGCGGAGCGCGACCCGAAAAGCTATGGTCTGTAACTTT"
    "TCGCGGGTCGAGCTAGTCCAAGTTCCGGCCTTTGTAATTCCGAAGTTGAATCGGTGATAC"
    "GGATTGACATGGGCCTAAACGTTCCGGCTGGTGTAGGATGATGCATCTCCAACATGTCTC"
    "TTACCGTTGCTGGGTCCGGCGGCTGTGGGATTGCGAGAGTGTCCGGCACCACCAATGTAC"
    "ACTTTCGGGAACACTCATTCGAAGAGGTTCTGCAGCTGCAGGCCTTGATACCTGCAGTCT"
    "GGGAGGCAATGCTGAGGCCCTCTGTTCCATGAAACCCGTACTATATCTTATGATGACAAT"
    "GAAATAGTCCTGTTTTACGACTCCAAGTTTCCTGCGCAATACCAAATACATTCCACGCGG"
    "CGCCTGGACTTAGTGTTCGT"
)

BASES = ["A", "C", "G", "T"]

# ─── Variant type distribution ─────────────────────────────────────────────────

# (type_name, n_train, n_test, ref_choice)
VARIANT_TYPES: list[tuple[str, int, int, str]] = [
    ("snv",           3500, 350, "snvindel"),
    ("indel_short",   2000, 200, "snvindel"),
    ("indel_medium",  1000, 100, "snvindel"),
    ("sv_del",         500,  50, "sv"),
    ("sv_dup",         300,  30, "sv"),
    ("sv_inv",         300,  30, "sv"),
    ("sv_large_del",   400,  40, "sv"),
    ("ins",           1000, 100, "sv"),
    ("invdel",         500,  50, "sv"),
]

TOTAL_TRAIN = sum(t[1] for t in VARIANT_TYPES)
TOTAL_TEST = sum(t[2] for t in VARIANT_TYPES)


# ─── Random DNA helpers ────────────────────────────────────────────────────────

def random_dna(rng: random.Random, length: int) -> str:
    """Return a random DNA string of the given length."""
    return "".join(rng.choice(BASES) for _ in range(length))


def ref_base_at(pos: int) -> str:
    """Return the reference base at a 1-based position in the snvindel ref.

    Falls back to 'N' if out of range.
    """
    idx = pos - 1
    if 0 <= idx < len(SNVINDEL_REF_SEQ):
        return SNVINDEL_REF_SEQ[idx]
    return "N"


def complement(base: str) -> str:
    """Return the complement of a DNA base."""
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(base, "N")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(complement(b) for b in reversed(seq))


# ─── Truth VCF generation ─────────────────────────────────────────────────────

def generate_snv_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate a single SNV variant record.

    Args:
        rng: Seeded random generator.
        pos: 1-based position on chr1.

    Returns:
        Dict with keys: chrom, pos, ref, alt, svtype_info.
    """
    ref = ref_base_at(pos)
    if ref == "N":
        ref = rng.choice(BASES)
    alt_choices = [b for b in BASES if b != ref]
    alt = rng.choice(alt_choices)
    return {"chrom": CHROM, "pos": pos, "ref": ref, "alt": alt, "info_extra": "TYPE=SNP"}


def generate_indel_short_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate a short indel (1-5bp insertion or deletion).

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    anchor = ref_base_at(pos)
    if anchor == "N":
        anchor = rng.choice(BASES)
    length = rng.randint(1, 5)
    if rng.random() < 0.5:
        # Insertion: REF=anchor, ALT=anchor+inserted_seq
        inserted = random_dna(rng, length)
        return {
            "chrom": CHROM, "pos": pos,
            "ref": anchor, "alt": anchor + inserted,
            "info_extra": f"TYPE=INDEL",
        }
    else:
        # Deletion: REF=anchor+deleted, ALT=anchor
        # Build ref from actual sequence if possible, else random
        deleted = ""
        for i in range(length):
            b = ref_base_at(pos + 1 + i)
            deleted += b if b != "N" else rng.choice(BASES)
        return {
            "chrom": CHROM, "pos": pos,
            "ref": anchor + deleted, "alt": anchor,
            "info_extra": "TYPE=INDEL",
        }


def generate_indel_medium_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate a medium indel (6-20bp insertion or deletion).

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    anchor = ref_base_at(pos)
    if anchor == "N":
        anchor = rng.choice(BASES)
    length = rng.randint(6, 20)
    if rng.random() < 0.5:
        inserted = random_dna(rng, length)
        return {
            "chrom": CHROM, "pos": pos,
            "ref": anchor, "alt": anchor + inserted,
            "info_extra": "TYPE=INDEL",
        }
    else:
        deleted = ""
        for i in range(length):
            b = ref_base_at(pos + 1 + i)
            deleted += b if b != "N" else rng.choice(BASES)
        return {
            "chrom": CHROM, "pos": pos,
            "ref": anchor + deleted, "alt": anchor,
            "info_extra": "TYPE=INDEL",
        }


def generate_sv_del_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate a small SV deletion (20-50bp).

    Uses 'N' as anchor REF following the SV VCF convention in this repo.

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(20, 50)
    deleted_seq = random_dna(rng, svlen)
    return {
        "chrom": CHROM, "pos": pos,
        "ref": "N" + deleted_seq, "alt": "N",
        "info_extra": f"SVTYPE=DEL;SVLEN={svlen}",
    }


def generate_sv_dup_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate an SV duplication (50-100bp).

    REF is the duplicated sequence; ALT is REF+REF (tandem duplication anchor).

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(50, 100)
    dup_seq = random_dna(rng, svlen)
    anchor = ref_base_at(pos)
    if anchor == "N":
        anchor = rng.choice(BASES)
    return {
        "chrom": CHROM, "pos": pos,
        "ref": anchor, "alt": anchor + dup_seq,
        "info_extra": f"SVTYPE=DUP;SVLEN={svlen}",
    }


def generate_sv_inv_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate an SV inversion (50-100bp).

    REF is the original sequence; ALT is the reverse complement.

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(50, 100)
    fwd_seq = random_dna(rng, svlen)
    rev_seq = reverse_complement(fwd_seq)
    return {
        "chrom": CHROM, "pos": pos,
        "ref": fwd_seq, "alt": rev_seq,
        "info_extra": f"SVTYPE=INV;SVLEN={svlen}",
    }


def generate_sv_large_del_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate a large SV deletion (100-200bp).

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(100, 200)
    deleted_seq = random_dna(rng, svlen)
    return {
        "chrom": CHROM, "pos": pos,
        "ref": "N" + deleted_seq, "alt": "N",
        "info_extra": f"SVTYPE=DEL;SVLEN={svlen}",
    }


def generate_ins_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate an SV insertion (70-100bp).

    Follows the same format as truth_ins_vaf*.vcf in this repo.

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(70, 100)
    anchor = ref_base_at(pos)
    if anchor == "N":
        anchor = rng.choice(BASES)
    inserted = random_dna(rng, svlen)
    return {
        "chrom": CHROM, "pos": pos,
        "ref": anchor, "alt": anchor + inserted,
        "info_extra": f"SVTYPE=INS;SVLEN={svlen}",
    }


def generate_invdel_variant(rng: random.Random, pos: int) -> dict[str, str]:
    """Generate an inversion-deletion variant (80-130bp).

    REF is the original sequence; ALT is the reverse complement of a shorter
    region (simulating deletion within the inverted segment), following the
    format in truth_invdel_vaf*.vcf.

    Args:
        rng: Seeded random generator.
        pos: 1-based anchor position.

    Returns:
        Dict with keys: chrom, pos, ref, alt, info_extra.
    """
    svlen = rng.randint(80, 130)
    # ALT is slightly shorter than REF (the deletion component)
    alt_len = svlen - rng.randint(5, 15)
    ref_seq = random_dna(rng, svlen)
    alt_seq = reverse_complement(random_dna(rng, alt_len))
    return {
        "chrom": CHROM, "pos": pos,
        "ref": ref_seq, "alt": alt_seq,
        "info_extra": f"SVTYPE=INVDEL;SVLEN={svlen}",
    }


def variant_generator_for(vtype: str):
    """Return the variant generator function for a given variant type.

    Args:
        vtype: Variant type string, e.g. 'snv', 'indel_short'.

    Returns:
        Callable(rng, pos) -> dict.
    """
    generators = {
        "snv": generate_snv_variant,
        "indel_short": generate_indel_short_variant,
        "indel_medium": generate_indel_medium_variant,
        "sv_del": generate_sv_del_variant,
        "sv_dup": generate_sv_dup_variant,
        "sv_inv": generate_sv_inv_variant,
        "sv_large_del": generate_sv_large_del_variant,
        "ins": generate_ins_variant,
        "invdel": generate_invdel_variant,
    }
    return generators[vtype]


def is_sv_type(vtype: str) -> bool:
    """Return True for variant types that require the SV reference."""
    return vtype in {"sv_del", "sv_dup", "sv_inv", "sv_large_del", "ins", "invdel"}


def ref_path_for(vtype: str) -> str:
    """Return the reference FASTA path (relative to repo root) for a type."""
    return SV_REF if is_sv_type(vtype) else SNVINDEL_REF


def write_truth_vcf(path: Path, vtype: str, variants: list[dict], vaf: float) -> None:
    """Write a VCF file containing the given variant records.

    Args:
        path: Output file path.
        vtype: Variant type (used to choose INFO fields).
        variants: List of variant dicts from the generator functions.
        vaf: True VAF for the INFO/VAF field.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={CHROM},length={CHROM_LEN}>\n")
        if is_sv_type(vtype):
            fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
            fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        else:
            fh.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">\n')
        fh.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            info = f"{v['info_extra']};VAF={vaf:.6f}"
            fh.write(
                f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t.\tPASS\t{info}\n"
            )


# ─── Parameter sampling ────────────────────────────────────────────────────────

def sample_params(idx: int) -> dict[str, Any]:
    """Sample varforge parameters for a given sample index.

    Uses a deterministic NumPy RNG seeded from the index so re-runs are
    reproducible and Nextflow caches remain valid.

    Args:
        idx: Global sample index (0-based, unique across train+test).

    Returns:
        Dict of sampled parameters.
    """
    rng = np.random.default_rng(idx * 7 + 42)
    coverage = float(np.exp(rng.uniform(np.log(1000), np.log(15000))))
    family_size_mean = float(rng.uniform(2.0, 8.0))
    pcr_cycles = int(rng.integers(5, 15))  # 5 inclusive, 15 exclusive
    fragment_mean = float(rng.uniform(130, 280))
    fragment_sd = float(rng.uniform(20, 50))
    mean_quality = int(rng.integers(30, 41))  # 30 inclusive, 41 exclusive
    vaf = float(np.exp(rng.uniform(np.log(0.0001), np.log(0.05))))
    seed = int(idx * 7 + 42)
    return {
        "coverage": round(coverage, 1),
        "family_size_mean": round(family_size_mean, 3),
        "pcr_cycles": pcr_cycles,
        "fragment_mean": round(fragment_mean, 1),
        "fragment_sd": round(fragment_sd, 1),
        "mean_quality": mean_quality,
        "vaf": vaf,
        "purity": min(vaf * 2, 1.0),  # diploid: purity = VAF * 2
        "seed": seed,
    }


# ─── Config YAML generation ────────────────────────────────────────────────────

def write_config(path: Path, sample_name: str, vtype: str, params: dict) -> None:
    """Write a varforge config YAML for a single sample.

    Args:
        path: Output YAML path.
        sample_name: Unique sample identifier.
        vtype: Variant type string.
        params: Parameter dict from sample_params().
    """
    split = "train" if "train" in str(path) else "test"
    ref = ref_path_for(vtype)
    sim_out = f"docs/benchmarking/ml-twist-duplex/simulations/{split}/{sample_name}"
    truth_vcf = f"docs/benchmarking/ml-twist-duplex/truth_vcfs/{sample_name}.vcf"

    yaml_text = f"""\
# varforge config: {sample_name}
reference: {ref}

output:
  directory: {sim_out}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {sample_name}
  read_length: 150
  coverage: {params['coverage']}
  platform: illumina

fragment:
  model: normal
  mean: {params['fragment_mean']}
  sd: {params['fragment_sd']}

quality:
  mean_quality: {params['mean_quality']}
  tail_decay: 0.002

tumour:
  purity: {params['purity']:.6f}
  ploidy: 2

mutations:
  vcf: {truth_vcf}

umi:
  length: 5
  duplex: true
  pcr_cycles: {params['pcr_cycles']}
  family_size_mean: {params['family_size_mean']}
  family_size_sd: 1.0
  inline: true

chromosomes:
  - {CHROM}

seed: {params['seed']}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml_text)


# ─── Sample name convention ────────────────────────────────────────────────────

def vaf_str(vaf: float) -> str:
    """Format VAF as a zero-padded percentage string for naming.

    Examples:
        >>> vaf_str(0.001)
        '0010pct'
        >>> vaf_str(0.05)
        '5000pct'
    """
    # Express as per-10000 to avoid floating-point noise
    pct_10000 = round(vaf * 10000)
    return f"{pct_10000:04d}pct"


def make_sample_name(vtype: str, vaf: float, idx: int, split: str) -> str:
    """Build a unique sample name.

    Format: {type}_{vaf_str}_{split_idx}
    Train indices: 0001, 0002, ...
    Test indices: t001, t002, ...

    Args:
        vtype: Variant type string.
        vaf: Sampled VAF.
        idx: Per-type index (1-based).
        split: 'train' or 'test'.

    Returns:
        Sample name string.
    """
    vs = vaf_str(vaf)
    if split == "train":
        return f"{vtype}_{vs}_{idx:04d}"
    return f"{vtype}_{vs}_t{idx:03d}"


# ─── Variant position selection ────────────────────────────────────────────────

def choose_positions(rng: random.Random, n_variants: int, vtype: str) -> list[int]:
    """Choose n_variants non-overlapping 1-based positions.

    Positions are spaced at least 150bp apart to avoid overlap in 150bp reads.
    For long alleles (SVs), a wider spacing is enforced.

    Args:
        rng: Seeded random generator.
        n_variants: Number of variants to place (1-3).
        vtype: Variant type, used to set minimum spacing.

    Returns:
        Sorted list of 1-based positions.
    """
    min_spacing = 300 if is_sv_type(vtype) else 150
    positions: list[int] = []
    attempts = 0
    max_attempts = 1000
    while len(positions) < n_variants and attempts < max_attempts:
        pos = rng.randint(MARGIN, MAX_POS)
        if all(abs(pos - p) >= min_spacing for p in positions):
            positions.append(pos)
        attempts += 1
    return sorted(positions)


# ─── Main generation logic ────────────────────────────────────────────────────

def build_sample_list() -> list[dict[str, Any]]:
    """Build the full list of 11,000 sample specifications.

    Returns:
        List of sample spec dicts, each containing name, vtype, split,
        params, config_path, vcf_path.
    """
    samples: list[dict[str, Any]] = []
    global_idx = 0

    for vtype, n_train, n_test, _ref_choice in VARIANT_TYPES:
        for split, n_samples in [("train", n_train), ("test", n_test)]:
            for local_idx in range(1, n_samples + 1):
                params = sample_params(global_idx)
                name = make_sample_name(vtype, params["vaf"], local_idx, split)

                if split == "train":
                    config_path = CONFIGS_TRAIN / f"{name}.yaml"
                else:
                    config_path = CONFIGS_TEST / f"{name}.yaml"

                vcf_path = TRUTH_VCFS_DIR / f"{name}.vcf"

                samples.append({
                    "name": name,
                    "vtype": vtype,
                    "split": split,
                    "global_idx": global_idx,
                    "local_idx": local_idx,
                    "params": params,
                    "config_path": config_path,
                    "vcf_path": vcf_path,
                })
                global_idx += 1

    return samples


def generate_sample(spec: dict[str, Any], dry_run: bool = False) -> None:
    """Generate the truth VCF and varforge config for a single sample.

    Args:
        spec: Sample specification dict from build_sample_list().
        dry_run: If True, skip writing files.
    """
    rng = random.Random(spec["global_idx"] * 13 + 7)
    vtype = spec["vtype"]
    params = spec["params"]
    vaf = params["vaf"]

    n_variants = rng.randint(1, 3)
    positions = choose_positions(rng, n_variants, vtype)
    generator = variant_generator_for(vtype)

    variants = []
    for pos in positions:
        try:
            v = generator(rng, pos)
            variants.append(v)
        except Exception as exc:
            print(
                f"[WARN] Variant generation failed for {spec['name']} at pos {pos}: {exc}",
                file=sys.stderr,
            )

    if not variants:
        print(f"[ERROR] No variants generated for {spec['name']}", file=sys.stderr)
        return

    if not dry_run:
        write_truth_vcf(spec["vcf_path"], vtype, variants, vaf)
        write_config(spec["config_path"], spec["name"], vtype, params)


def write_manifest(samples: list[dict[str, Any]], dry_run: bool = False) -> None:
    """Write the dataset manifest JSON.

    The manifest records every sample's split, variant type, and parameters
    so downstream scripts can cross-reference without re-reading YAML files.

    Args:
        samples: List of sample spec dicts.
        dry_run: If True, skip writing.
    """
    manifest: dict[str, Any] = {
        "total": len(samples),
        "n_train": sum(1 for s in samples if s["split"] == "train"),
        "n_test": sum(1 for s in samples if s["split"] == "test"),
        "samples": [],
    }

    for s in samples:
        p = s["params"]
        manifest["samples"].append({
            "name": s["name"],
            "split": s["split"],
            "vtype": s["vtype"],
            "global_idx": s["global_idx"],
            "coverage": p["coverage"],
            "family_size_mean": p["family_size_mean"],
            "pcr_cycles": p["pcr_cycles"],
            "fragment_mean": p["fragment_mean"],
            "fragment_sd": p["fragment_sd"],
            "mean_quality": p["mean_quality"],
            "vaf": p["vaf"],
            "purity": p["purity"],
            "seed": p["seed"],
            "config": str(s["config_path"].relative_to(REPO)),
            "truth_vcf": str(s["vcf_path"].relative_to(REPO)),
        })

    if not dry_run:
        manifest_path = ML_DIR / "manifest.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        print(f"Manifest written: {manifest_path}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate varforge configs and truth VCFs for the Twist duplex ML dataset.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be generated without writing files.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Generate only the first N samples (for testing).",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    print(f"Building sample list ({TOTAL_TRAIN} train + {TOTAL_TEST} test = {TOTAL_TRAIN + TOTAL_TEST} total)...")
    samples = build_sample_list()

    if args.limit is not None:
        samples = samples[: args.limit]
        print(f"[DRY-RUN] Limiting to {len(samples)} samples.")

    if args.dry_run:
        for s in samples[:5]:
            print(f"  {s['split']:5s}  {s['vtype']:15s}  {s['name']}")
        print(f"  ... ({len(samples)} total)")
        return

    total = len(samples)
    for i, spec in enumerate(samples):
        generate_sample(spec, dry_run=args.dry_run)
        if (i + 1) % 500 == 0 or (i + 1) == total:
            print(f"  [{i + 1}/{total}] generated", flush=True)

    write_manifest(samples, dry_run=args.dry_run)

    # Print summary
    type_counts: dict[str, dict[str, int]] = {}
    for s in samples:
        type_counts.setdefault(s["vtype"], {"train": 0, "test": 0})
        type_counts[s["vtype"]][s["split"]] += 1

    print(f"\nGenerated {total} samples:")
    print(f"  {'Type':<16} {'Train':>6} {'Test':>6}")
    for vtype, counts in type_counts.items():
        print(f"  {vtype:<16} {counts['train']:>6} {counts['test']:>6}")

    n_vcfs = sum(1 for s in samples if s["vcf_path"].exists())
    n_configs = sum(
        1 for s in samples
        if s["config_path"].exists()
    )
    print(f"\nTruth VCFs written:    {n_vcfs}")
    print(f"Config YAMLs written:  {n_configs}")


if __name__ == "__main__":
    main()
