#!/usr/bin/env python3
"""Generate varforge configs and truth VCFs for fusion/translocation benchmarking.

Produces 25 VAF levels × 2 replicates = 50 configs and 50 truth VCFs,
all targeting the sv_new benchmark suite.

## Approach

varforge does not natively support BND/translocation variants in its VCF input
parser (BND records are explicitly skipped). The alternative used here is the
"fused reference" approach:

  1. A synthetic reference (fusion_ref.fa) is constructed by concatenating
     the gene A segment (chr1:0-200) with the gene B segment (chr1:900-2000).
     The breakpoint falls at position 200 in this fused reference.

  2. varforge simulates reads from fusion_ref.fa without any mutations VCF.
     These reads represent fusion-bearing molecules: reads crossing position 200
     span the breakpoint junction and contain k-mers unique to the fusion.

  3. At simulation time, purity controls the fraction of fusion-bearing reads.
     Since each simulated molecule in this reference is a fusion molecule,
     purity = 2 × VAF (diploid, one fusion allele out of two chromosomes).

  4. The resulting FASTQ files are mixed with wild-type reads from the normal
     chr1 reference to produce a final mixed sample. The mixing is done by the
     downstream benchmark runner (see README.md for the mixing step).

The fusion target for kam is in fusion_targets.fa:
  BCR_ABL1__chr1:150-200__chr1:900-950__fusion
  (50bp of gene A ending at the breakpoint, 50bp of gene B starting there)

The truth VCF uses BND notation. For benchmarking purposes, a call is a true
positive if kam detects the BCR_ABL1 fusion at the correct breakpoint.

Run from the repo root:
    python3 docs/benchmarking/sv_new/scripts/make_fusion_suite.py
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

# Fused reference: gene A (chr1:0-200) concatenated with gene B (chr1:900-2000).
# Breakpoint is at position 200 in this reference.
FUSION_REF_PATH = "docs/benchmarking/sv_new/data/fusion_ref.fa"

# Normal (wild-type) reference for depth normalisation.
NORMAL_REF_PATH = "docs/benchmarking/sv/data/ref.fa"

# Fusion breakpoint coordinates in the original chr1 coordinate space.
# Gene A: chr1 position 200 (0-based), gene B: chr1 position 900 (0-based).
# VCF POS fields are 1-based.
BREAKPOINT_A_POS = 200  # 1-based VCF position
BREAKPOINT_B_POS = 900  # 1-based VCF position

# Anchor bases at each breakpoint (from the chr1 reference sequence).
# ref.fa seq[199] = 'A' (1-based position 200), seq[899] = 'T' (1-based position 900)
ANCHOR_BASE_A = "A"
ANCHOR_BASE_B = "T"

VCF_HEADERS = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=2000>",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate BND record ID">',
    '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]

# varforge config for the fusion reference. No mutations VCF: all simulated
# reads are fusion-bearing by design (they come from the fused reference).
# Purity = 2 × VAF (diploid; one fusion allele in a two-allele background).
# The downstream mixing step blends these reads with wild-type reads.
CONFIG_TEMPLATE = """\
# varforge config: {sample_name} — VAF {vaf_pct:.3f}%
# Fusion simulation using the fused reference approach.
# All reads from this simulation cross the BCR_ABL1 breakpoint (chr1:200|900).
# Mix with wild-type reads (from the normal ref) at ratio purity:{wt_ratio} to achieve
# the target fusion VAF. See docs/benchmarking/sv_new/README.md for mixing steps.
reference: {ref_path}

output:
  directory: {out_dir}
  fastq: true
  bam: false
  truth_vcf: true
  manifest: true

sample:
  name: {sample_name}
  read_length: 150
  coverage: {fusion_coverage:.1f}
  platform: illumina

fragment:
  model: normal
  mean: 167.0
  sd: 30.0

quality:
  mean_quality: 37
  tail_decay: 0.002

tumour:
  purity: 1.0
  ploidy: 2

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

# Wild-type background config (normal chr1, no mutations).
WT_CONFIG_TEMPLATE = """\
# varforge config: {sample_name}_WT — wild-type background for fusion mixing.
reference: {ref_path}

output:
  directory: {out_dir}_wt
  fastq: true
  bam: false
  truth_vcf: false
  manifest: false

sample:
  name: {sample_name}_WT
  read_length: 150
  coverage: {wt_coverage:.1f}
  platform: illumina

fragment:
  model: normal
  mean: 167.0
  sd: 30.0

quality:
  mean_quality: 37
  tail_decay: 0.002

tumour:
  purity: 0.0
  ploidy: 2

umi:
  length: 5
  duplex: true
  pcr_cycles: 8
  family_size_mean: 4.0
  family_size_sd: 1.5
  inline: true

chromosomes:
  - chr1

seed: {wt_seed}
"""


def vaf_tag(vaf: float) -> str:
    """Return zero-padded 4-digit VAF tag (e.g. 0.01 → '0100')."""
    return f"{round(vaf * 10000):04d}"


def write_truth_vcf(path: Path, vaf: float) -> None:
    """Write a BND truth VCF for one fusion sample at the given VAF.

    The fusion joins chr1:200 (gene A breakpoint) to chr1:900 (gene B
    breakpoint). Two BND records are emitted: one anchored at each breakpoint.
    """
    fusion_id = "BCR_ABL1"
    # Record 1: gene A side — forward-forward orientation.
    rec1 = (
        f"chr1\t{BREAKPOINT_A_POS}\t{fusion_id}_1\t{ANCHOR_BASE_A}"
        f"\t{ANCHOR_BASE_A}]chr1:{BREAKPOINT_B_POS}]"
        f"\t.\tPASS"
        f"\tSVTYPE=BND;MATEID={fusion_id}_2;VAF={vaf}"
    )
    # Record 2: gene B side — forward-forward orientation.
    rec2 = (
        f"chr1\t{BREAKPOINT_B_POS}\t{fusion_id}_2\t{ANCHOR_BASE_B}"
        f"\t]chr1:{BREAKPOINT_A_POS}]{ANCHOR_BASE_B}"
        f"\t.\tPASS"
        f"\tSVTYPE=BND;MATEID={fusion_id}_1;VAF={vaf}"
    )
    with open(path, "w") as f:
        for h in VCF_HEADERS:
            f.write(h + "\n")
        f.write(rec1 + "\n")
        f.write(rec2 + "\n")


def write_fusion_config(
    path: Path,
    sample_name: str,
    vaf: float,
    out_dir: str,
    seed: int,
) -> None:
    """Write a varforge YAML config for one fusion replicate.

    Total coverage is 5000×. At VAF=vaf, fusion reads make up vaf of the total.
    The fusion simulation generates reads at coverage = 5000 * vaf from the
    fused reference. The wild-type config generates the remaining background.
    """
    total_coverage = 5000.0
    fusion_coverage = total_coverage * vaf
    wt_coverage = total_coverage * (1.0 - vaf)
    # Express the mixing ratio as integers for the comment.
    pct_int = round(vaf * 100, 4)
    wt_ratio_int = round((1.0 - vaf) * 100, 4)

    with open(path, "w") as f:
        f.write(CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            vaf_pct=vaf * 100,
            ref_path=FUSION_REF_PATH,
            out_dir=out_dir,
            fusion_coverage=fusion_coverage,
            wt_ratio=wt_ratio_int,
            seed=seed,
        ))


def write_wt_config(
    path: Path,
    sample_name: str,
    vaf: float,
    out_dir: str,
    seed: int,
) -> None:
    """Write a varforge YAML config for the wild-type background reads."""
    total_coverage = 5000.0
    wt_coverage = total_coverage * (1.0 - vaf)

    with open(path, "w") as f:
        f.write(WT_CONFIG_TEMPLATE.format(
            sample_name=sample_name,
            ref_path=NORMAL_REF_PATH,
            out_dir=out_dir,
            wt_coverage=wt_coverage,
            wt_seed=seed + 500,
        ))


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CFGS.mkdir(parents=True, exist_ok=True)

    count = 0
    for vaf in VAF_LEVELS:
        t = vaf_tag(vaf)
        tag_int = round(vaf * 10000)

        for rep_idx, rep in enumerate(["a", "b"]):
            rep_up = rep.upper()
            seed_offset = rep_idx * 1000

            sample_name = f"FUSION_VAF{t}_{rep_up}"
            vcf_path = DATA / f"truth_fusion_vaf{t}_{rep}.vcf"
            cfg_path = CFGS / f"fusion_vaf{t}_{rep}.yaml"
            wt_cfg_path = CFGS / f"fusion_vaf{t}_{rep}_wt.yaml"
            out_dir = f"docs/benchmarking/sv_new/results/sim_fusion_vaf{t}_{rep}"
            # Seed range 130000+ to avoid collisions with existing suites.
            seed = 130000 + tag_int + seed_offset

            write_truth_vcf(vcf_path, vaf)
            write_fusion_config(cfg_path, sample_name, vaf, out_dir, seed)
            write_wt_config(wt_cfg_path, sample_name, vaf, out_dir, seed)
            count += 1

    print(f"Generated {count} fusion config pairs and truth VCFs in {ROOT}/")
    print()
    print("Each VAF level produces three files:")
    print("  configs/fusion_vaf<T>_<rep>.yaml      — fusion reads (from fused reference)")
    print("  configs/fusion_vaf<T>_<rep>_wt.yaml   — wild-type background reads")
    print("  data/truth_fusion_vaf<T>_<rep>.vcf    — BND truth VCF")
    print()
    print("See docs/benchmarking/sv_new/README.md for the mixing step.")


if __name__ == "__main__":
    main()
