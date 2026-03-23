#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for the SNV/indel benchmark suite.

Produces 25 VAF levels × 2 replicates × 2 types (SNV, indel) = 100 configs +
100 truth VCFs. Run from the repo root.
"""
import re
from pathlib import Path

# VAF levels as fractions (0.05% = 0.0005)
VAF_LEVELS = [
    0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0050,
    0.0060, 0.0075, 0.0100, 0.0125, 0.0150, 0.0175, 0.0200, 0.0250, 0.0300,
    0.0350, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.1000,
]

ROOT = Path("docs/benchmarking/snvindel")
DATA  = ROOT / "data"
CFGS  = ROOT / "configs"

SNV_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=2000>
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

INDEL_HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=2000>
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Indel length (negative=deletion)">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

SNV_ROWS = [
    ("chr1", 50,   "A", "T",      "TYPE=SNP"),
    ("chr1", 200,  "A", "G",      "TYPE=SNP"),
    ("chr1", 450,  "G", "A",      "TYPE=SNP"),
    ("chr1", 800,  "G", "A",      "TYPE=SNP"),
    ("chr1", 1200, "A", "C",      "TYPE=SNP"),
]

INDEL_ROWS = [
    ("chr1", 100,  "TTTG",  "T",      "TYPE=INDEL;SVLEN=-3"),
    ("chr1", 300,  "A",     "AGACGT", "TYPE=INDEL;SVLEN=5"),
    ("chr1", 550,  "GTGA",  "G",      "TYPE=INDEL;SVLEN=-3"),
    ("chr1", 750,  "T",     "TGCTAG", "TYPE=INDEL;SVLEN=5"),
    ("chr1", 1000, "TACA",  "T",      "TYPE=INDEL;SVLEN=-3"),
]


def vaf_tag(vaf: float) -> str:
    return f"{round(vaf * 10000):04d}"


def write_vcf(path: Path, header: str, rows: list, vaf: float) -> None:
    with open(path, "w") as f:
        f.write(header)
        for chrom, pos, ref, alt, info in rows:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info};VAF={vaf}\n")


CONFIG_TEMPLATE = """\
# varforge config: {sample_name} — VAF {vaf_pct:.3f}%
reference: docs/benchmarking/snvindel/data/ref.fa

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


def write_config(path: Path, sample_name: str, vaf: float, purity: float,
                 vcf_rel: str, out_dir: str, seed: int) -> None:
    with open(path, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            vaf_pct=vaf * 100,
            purity=purity,
            vcf_path=vcf_rel,
            out_dir=out_dir,
            seed=seed,
        ))


n_snv = n_indel = 0
for vaf in VAF_LEVELS:
    t = vaf_tag(vaf)
    purity = round(vaf * 2, 6)
    tag_int = round(vaf * 10000)

    for rep_idx, rep in enumerate(["a", "b"]):
        rep_up = rep.upper()
        seed_offset = rep_idx * 1000

        # --- SNV ---
        snv_vcf  = DATA / f"truth_snvs_vaf{t}_{rep}.vcf"
        snv_cfg  = CFGS / f"snv_vaf{t}_{rep}.yaml"
        snv_seed = 10000 + tag_int + seed_offset
        write_vcf(snv_vcf, SNV_HEADER, SNV_ROWS, vaf)
        write_config(
            snv_cfg,
            sample_name=f"SNV_VAF{t}_{rep_up}",
            vaf=vaf, purity=purity,
            vcf_rel=f"docs/benchmarking/snvindel/data/truth_snvs_vaf{t}_{rep}.vcf",
            out_dir=f"docs/benchmarking/snvindel/results/sim_snv_vaf{t}_{rep}",
            seed=snv_seed,
        )
        n_snv += 1

        # --- Indel ---
        indel_vcf  = DATA / f"truth_indels_vaf{t}_{rep}.vcf"
        indel_cfg  = CFGS / f"indel_vaf{t}_{rep}.yaml"
        indel_seed = 30000 + tag_int + seed_offset
        write_vcf(indel_vcf, INDEL_HEADER, INDEL_ROWS, vaf)
        write_config(
            indel_cfg,
            sample_name=f"INDEL_VAF{t}_{rep_up}",
            vaf=vaf, purity=purity,
            vcf_rel=f"docs/benchmarking/snvindel/data/truth_indels_vaf{t}_{rep}.vcf",
            out_dir=f"docs/benchmarking/snvindel/results/sim_indel_vaf{t}_{rep}",
            seed=indel_seed,
        )
        n_indel += 1

print(f"Generated {n_snv} SNV configs + {n_indel} indel configs ({n_snv + n_indel} total)")
