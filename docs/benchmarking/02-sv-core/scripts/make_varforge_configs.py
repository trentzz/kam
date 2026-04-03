#!/usr/bin/env python3
"""Generate varforge simulation configs for SV benchmarking.

Creates one config per VAF level (0.5%, 1%, 2%, 5%).  Each config:
  - Uses the synthetic 2000 bp reference.
  - Injects SVs via VCF (full-sequence alleles for DEL/DUP, and INV).
  - Uses Twist duplex UMI chemistry: length 5, inline, duplex, skip 2.
  - cfDNA fragment model: mean 167 bp, SD 30 bp.
  - Read length 150 bp.
  - Coverage 5000x raw (≈1000 molecule families) to stress-test low VAF.
  - Restricted to chr1.
"""

import os

BASE_DIR = "docs/benchmarking/sv"
DATA_DIR = f"{BASE_DIR}/data"
CONFIGS_DIR = f"{BASE_DIR}/configs"
REF_PATH = f"{DATA_DIR}/ref.fa"

VAFS = [
    (0.005, "005"),
    (0.010, "010"),
    (0.020, "020"),
    (0.050, "050"),
]

TEMPLATE = """\
# varforge config: SV benchmarking at VAF {vaf_pct}%
reference: {ref}

output:
  directory: {out_dir}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: SV_VAF{vaf_tag}
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
  purity: {purity}
  ploidy: 2

mutations:
  vcf: {vcf}

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


def main() -> None:
    os.makedirs(CONFIGS_DIR, exist_ok=True)
    for vaf, tag in VAFS:
        # Convert VAF to purity: purity * 0.5 = VAF for diploid (het variant).
        # Simpler: use purity=1.0 and set VAF directly in truth VCF INFO.
        # varforge uses purity to scale VAF.  Since INFO VAF overrides purity
        # in custom_mutations mode, set purity=1.0 and let the VCF VAF field drive.
        purity = vaf * 2  # diploid het: effective VAF = purity/2
        vcf_path = f"{DATA_DIR}/truth_svs_vaf{tag}.vcf"
        out_dir = f"{BASE_DIR}/results/sim_vaf{tag}"
        cfg_content = TEMPLATE.format(
            vaf_pct=vaf * 100,
            vaf_tag=tag,
            ref=REF_PATH,
            out_dir=out_dir,
            purity=purity,
            vcf=vcf_path,
            seed=1000 + int(tag),
        )
        cfg_path = f"{CONFIGS_DIR}/sim_vaf{tag}.yaml"
        with open(cfg_path, "w") as f:
            f.write(cfg_content)
        print(f"Wrote config: {cfg_path}")


if __name__ == "__main__":
    main()
