#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for InvDel benchmarking.

Produces 25 VAF levels × 2 replicates = 50 configs and 50 truth VCFs,
all targeting the sv_new benchmark suite.

The InvDel variant is an 80 bp event fully contained within the 200 bp
target window chr1:801-1000 on the 2000 bp synthetic reference:
  - REF: ref[860:940] (1-based POS=861), 80 bp
  - ALT: RC(ref[860:920]), 60 bp (inversion of first 60 bp; last 20 bp deleted)
  - Target window: ref[800:1000] (1-based 801-1000), 200 bp
  - Junction: last 35 bp of ALT + first 35 bp of ref after variant (ref[940:975])

The variant is fully inside the target window (0-based [860,940) inside [800,1000)),
so no part of the REF allele extends before or after the window.

Run from the repo root:
    python3 docs/benchmarking/03-sv-extended/scripts/make_invdel_suite.py
"""

from pathlib import Path

VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

ROOT = Path("docs/benchmarking/sv_new")
DATA = ROOT / "data"
CFGS = ROOT / "configs"

# The reference is shared with all other SV benchmarks.
REF_PATH = "docs/benchmarking/sv/data/ref.fa"

# InvDel variant definition: 80 bp event at chr1:861 (1-based).
# REF: ref[860:940] (0-based), 80 bp, fully within target window ref[800:1000].
# ALT: RC(ref[860:920]), 60 bp. First 60 bp are inverted; last 20 bp are deleted.
# alt path (60 bp) < ref path (80 bp) confirms the deletion component.
# Verified: REF matches ref.fa exactly; ALT is RC of the inverted segment.
INVDEL_POS = 861
INVDEL_HEADERS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=2000>",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
    '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]

# REF: ref[860:940] from the 2000 bp synthetic reference (80 bp).
INVDEL_REF = "GAGGGTCCTCCCATCTCCTGTGATGCATGGTGTGCTTACTGGGATGAATGCGCCGCAAGTAGCAGGTCCCGGCGTGGATA"
# ALT: RC(ref[860:920]) — inversion of the first 60 bp, with the last 20 bp deleted (60 bp).
INVDEL_ALT = "ACTTGCGGCGCATTCATCCCAGTAAGCACACCATGCATCACAGGAGATGGGAGGACCCTC"


CONFIG_TEMPLATE = """\
# varforge config: {sample_name} — VAF {vaf_pct:.3f}%
reference: docs/benchmarking/sv/data/ref.fa

output:
  directory: {out_dir}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {sample_name}
  read_length: 150
  coverage: 5000.0
  platform: illumina

fragment:
  model: normal
  mean: 167.0
  sd: 30.0

quality:
  mean_quality: 37
  tail_decay: 0.002

tumour:
  purity: {purity:.6f}
  ploidy: 2

mutations:
  vcf: {vcf_path}

umi:
  length: 5
  duplex: true
  pcr_cycles: 8
  family_size_mean: 4.0
  family_size_sd: 1.5
  inline: true

chromosomes:
  - chr1

seed: {seed}
"""


def vaf_tag(vaf: float) -> str:
    """Return zero-padded 4-digit VAF tag (e.g. 0.01 → '0100')."""
    return f"{round(vaf * 10000):04d}"


def write_vcf(path: Path, vaf: float) -> None:
    """Write a truth VCF for one InvDel sample at the given VAF."""
    row = (
        f"chr1\t{INVDEL_POS}\t.\t{INVDEL_REF}\t{INVDEL_ALT}\t.\tPASS"
        f"\tSVTYPE=INVDEL;SVLEN=80;VAF={vaf}"
    )
    with open(path, "w") as f:
        for h in INVDEL_HEADERS:
            f.write(h + "\n")
        f.write(row + "\n")


def write_config(
    path: Path,
    sample_name: str,
    vaf: float,
    purity: float,
    vcf_rel: str,
    out_dir: str,
    seed: int,
) -> None:
    """Write a varforge YAML config for one InvDel sample."""
    with open(path, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            vaf_pct=vaf * 100,
            purity=purity,
            vcf_path=vcf_rel,
            out_dir=out_dir,
            seed=seed,
        ))


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CFGS.mkdir(parents=True, exist_ok=True)

    count = 0
    for vaf in VAF_LEVELS:
        t = vaf_tag(vaf)
        purity = round(vaf * 2, 6)
        tag_int = round(vaf * 10000)

        for rep_idx, rep in enumerate(["a", "b"]):
            rep_up = rep.upper()
            seed_offset = rep_idx * 1000

            sample_name = f"INVDEL_VAF{t}_{rep_up}"
            vcf_path = DATA / f"truth_invdel_vaf{t}_{rep}.vcf"
            cfg_path = CFGS / f"invdel_vaf{t}_{rep}.yaml"
            vcf_rel = f"docs/benchmarking/03-sv-extended/data/truth_invdel_vaf{t}_{rep}.vcf"
            out_dir = f"docs/benchmarking/03-sv-extended/results/sim_invdel_vaf{t}_{rep}"
            # Seed range 90000+ mirrors the existing sv/ suite convention.
            seed = 90000 + tag_int + seed_offset

            write_vcf(vcf_path, vaf)
            write_config(cfg_path, sample_name, vaf, purity, vcf_rel, out_dir, seed)
            count += 1

    print(f"Generated {count} InvDel configs and truth VCFs in {ROOT}/")


if __name__ == "__main__":
    main()
