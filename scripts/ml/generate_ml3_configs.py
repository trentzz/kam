"""Generate varforge configs for the ML3 dataset with an explicit train/test split.

Train split: 10,000 configs
Test split:  1,000 configs

Configs are written to separate subdirectories so the pipeline runner can
process them independently and prevent any data leakage between splits.

Output layout:
  docs/benchmarking/ml/configs/train/   — 10k YAML + params.json pairs
  docs/benchmarking/ml/configs/test/    — 1k YAML + params.json pairs
  docs/benchmarking/ml/configs/ml3_train_manifest.json
  docs/benchmarking/ml/configs/ml3_test_manifest.json

Each config has a companion <name>_params.json recording the simulation
parameters so the training script can use them as features.
"""

import json
import math
import random
from pathlib import Path

REPO = Path("/home/trent/code/kam")
CONFIGS_DIR = REPO / "docs/benchmarking/ml/configs"
TRAIN_DIR = CONFIGS_DIR / "train"
TEST_DIR = CONFIGS_DIR / "test"

# Available VAF tags in the truth VCF filenames (integer = vaf * 10000).
AVAILABLE_VAF_TAGS = [
    5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75,
    100, 125, 150, 175, 200, 250, 300, 350, 400,
    500, 600, 700, 800, 1000,
]

# Reference files per variant type.
REFERENCES = {
    "snv":    "docs/benchmarking/snvindel/data/ref.fa",
    "indel":  "docs/benchmarking/snvindel/data/ref.fa",
    "sv":     "docs/benchmarking/sv/data/ref.fa",
    "ins":    "docs/benchmarking/sv/data/ref.fa",
    "invdel": "docs/benchmarking/sv/data/ref.fa",
}

# Truth VCF path templates — {vaf} = 4-digit tag, {rep} = a or b.
TRUTH_VCF_TEMPLATES = {
    "snv":    "docs/benchmarking/snvindel/data/truth_snvs_vaf{vaf}_{rep}.vcf",
    "indel":  "docs/benchmarking/snvindel/data/truth_indels_vaf{vaf}_{rep}.vcf",
    "sv":     "docs/benchmarking/sv/data/truth_svs_vaf{vaf}_{rep}.vcf",
    "ins":    "docs/benchmarking/sv/data/truth_ins_vaf{vaf}_{rep}.vcf",
    "invdel": "docs/benchmarking/sv/data/truth_invdel_vaf{vaf}_{rep}.vcf",
}

# Per-type counts: balanced across all five types.
TRAIN_COUNTS = {
    "snv":    2000,
    "indel":  2000,
    "sv":     2000,
    "ins":    2000,
    "invdel": 2000,
}

TEST_COUNTS = {
    "snv":    200,
    "indel":  200,
    "sv":     200,
    "ins":    200,
    "invdel": 200,
}


def nearest_vaf_tag(vaf_fraction: float) -> int:
    """Return the closest available VAF tag (integer = vaf * 10000)."""
    vaf_int = vaf_fraction * 10000.0
    return min(AVAILABLE_VAF_TAGS, key=lambda t: abs(t - vaf_int))


def sample_params(rng: random.Random) -> dict:
    """Sample simulation parameters with log-uniform coverage and VAF."""
    coverage = 10 ** rng.uniform(math.log10(500), math.log10(20000))
    vaf = 10 ** rng.uniform(math.log10(0.0003), math.log10(0.15))
    family_size_mean = 10 ** rng.uniform(math.log10(1.5), math.log10(10.0))
    family_size_sd = family_size_mean * rng.uniform(0.2, 0.6)
    pcr_cycles = rng.randint(4, 14)
    fragment_mean = rng.uniform(130, 280)
    fragment_sd = rng.uniform(15, 50)
    mean_quality = rng.uniform(30, 40)
    tail_decay = rng.uniform(0.001, 0.005)
    purity = min(vaf * 2.0, 1.0)
    return {
        "coverage": coverage,
        "vaf": vaf,
        "family_size_mean": family_size_mean,
        "family_size_sd": family_size_sd,
        "pcr_cycles": pcr_cycles,
        "fragment_mean": fragment_mean,
        "fragment_sd": fragment_sd,
        "mean_quality": mean_quality,
        "tail_decay": tail_decay,
        "purity": purity,
        "ploidy": 2,
    }


def config_seed(name: str) -> int:
    """Deterministic seed from config name. Offset from ml2 range (200000+)."""
    return hash(name) % 100000 + 400000


def make_config_yaml(
    name: str,
    vtype: str,
    params: dict,
    truth_vcf: str,
    output_dir: str,
) -> str:
    """Produce the varforge YAML config string."""
    vaf_tag = nearest_vaf_tag(params["vaf"])
    vaf_label = f"{vaf_tag / 10000:.4%}"
    seed = config_seed(name)
    ref = REFERENCES[vtype]

    lines = [
        f"# varforge config: {name.upper()} — VAF {vaf_label}",
        f"reference: {ref}",
        "",
        "output:",
        f"  directory: {output_dir}",
        "  fastq: true",
        "  bam: false",
        "  truth_vcf: true",
        "  manifest: true",
        "",
        "sample:",
        f"  name: {name.upper()}",
        "  read_length: 150",
        f"  coverage: {params['coverage']:.1f}",
        "  platform: illumina",
        "",
        "fragment:",
        "  model: normal",
        f"  mean: {params['fragment_mean']:.1f}",
        f"  sd: {params['fragment_sd']:.1f}",
        "",
        "quality:",
        f"  mean_quality: {int(round(params['mean_quality']))}",
        f"  tail_decay: {params['tail_decay']:.4f}",
        "",
        "tumour:",
        f"  purity: {params['purity']:.6f}",
        f"  ploidy: {params['ploidy']}",
        "",
        "mutations:",
        f"  vcf: {truth_vcf}",
        "",
        "umi:",
        "  length: 5",
        "  duplex: true",
        f"  pcr_cycles: {params['pcr_cycles']}",
        f"  family_size_mean: {params['family_size_mean']:.2f}",
        f"  family_size_sd: {params['family_size_sd']:.2f}",
        "",
        "chromosomes:",
        "  - chr1",
        "",
        f"seed: {seed}",
    ]
    return "\n".join(lines) + "\n"


def make_params_json(
    name: str,
    vtype: str,
    params: dict,
    replicate: str,
    split: str,
) -> dict:
    """Produce the params.json dict stored alongside each config."""
    vaf_tag = nearest_vaf_tag(params["vaf"])
    return {
        "coverage": round(params["coverage"], 2),
        "family_size_mean": round(params["family_size_mean"], 3),
        "pcr_cycles": params["pcr_cycles"],
        "fragment_mean": round(params["fragment_mean"], 2),
        "variant_type": vtype,
        "vaf": round(params["vaf"], 6),
        "param_variant": "ml3",
        "replicate": replicate,
        "vaf_tag": vaf_tag,
        "split": split,
    }


def generate_split(
    split: str,
    type_counts: dict,
    out_dir: Path,
    rng: random.Random,
) -> list[str]:
    """Generate configs for one split (train or test). Returns list of names."""
    out_dir.mkdir(parents=True, exist_ok=True)
    names = []

    for vtype, count in type_counts.items():
        for idx in range(1, count + 1):
            p = sample_params(rng)
            replicate = rng.choice(["a", "b"])
            vaf_tag = nearest_vaf_tag(p["vaf"])
            vaf_tag_str = f"{vaf_tag:04d}"
            vaf_label = int(round(p["vaf"] * 10000))

            name = f"{vtype}_ml3_{split}_vaf{vaf_label:04d}_{idx:04d}"

            truth_vcf = TRUTH_VCF_TEMPLATES[vtype].format(
                vaf=vaf_tag_str, rep=replicate
            )
            output_dir = f"docs/benchmarking/ml/results/{split}/sim_{name}"

            yaml_content = make_config_yaml(name, vtype, p, truth_vcf, output_dir)
            (out_dir / f"{name}.yaml").write_text(yaml_content)

            params_data = make_params_json(name, vtype, p, replicate, split)
            (out_dir / f"{name}_params.json").write_text(
                json.dumps(params_data, indent=2) + "\n"
            )

            names.append(name)

    return names


def main() -> None:
    """Generate all train and test configs."""
    # Use a fixed top-level seed; train and test use non-overlapping sub-sequences.
    rng_train = random.Random(1001)
    rng_test = random.Random(1002)

    print("Generating training configs...")
    train_names = generate_split("train", TRAIN_COUNTS, TRAIN_DIR, rng_train)
    manifest_train = CONFIGS_DIR / "ml3_train_manifest.json"
    manifest_train.write_text(json.dumps(train_names, indent=2) + "\n")
    print(f"  {len(train_names)} training configs written to {TRAIN_DIR}")

    print("Generating test configs...")
    test_names = generate_split("test", TEST_COUNTS, TEST_DIR, rng_test)
    manifest_test = CONFIGS_DIR / "ml3_test_manifest.json"
    manifest_test.write_text(json.dumps(test_names, indent=2) + "\n")
    print(f"  {len(test_names)} test configs written to {TEST_DIR}")

    print(f"\nTotal: {len(train_names)} train + {len(test_names)} test configs.")


if __name__ == "__main__":
    main()
