"""Generate 3000 diverse varforge configs for the kam ML training dataset (generation 2).

Uses random sampling across the full parameter space rather than grid sweeps,
producing a richer and more diverse training set for better generalisation.

Output:
  bigdata/experiments/02-ml-single-strand/configs/  — one YAML + one params.json per config
  bigdata/experiments/02-ml-single-strand/configs/ml2_manifest.json  — list of all 3000 names
"""

import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONFIGS_DIR = REPO / "bigdata/experiments/02-ml-single-strand/configs"

# Available VAF tags in the truth VCF filenames (as integers, vaf * 10000).
AVAILABLE_VAF_TAGS = [
    5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75,
    100, 125, 150, 175, 200, 250, 300, 350, 400,
    500, 600, 700, 800, 1000,
]

# Reference files per variant type.
REFERENCES = {
    "snv":    "docs/benchmarking/01-snvindel/data/ref.fa",
    "indel":  "docs/benchmarking/01-snvindel/data/ref.fa",
    "sv":     "docs/benchmarking/02-sv-core/data/ref.fa",
    "ins":    "docs/benchmarking/02-sv-core/data/ref.fa",
    "invdel": "docs/benchmarking/02-sv-core/data/ref.fa",
}

# Truth VCF path templates per variant type.
# {vaf} will be replaced with the 4-digit VAF tag, {rep} with a or b.
TRUTH_VCF_TEMPLATES = {
    "snv":    "docs/benchmarking/01-snvindel/data/truth_snvs_vaf{vaf}_{rep}.vcf",
    "indel":  "docs/benchmarking/01-snvindel/data/truth_indels_vaf{vaf}_{rep}.vcf",
    "sv":     "docs/benchmarking/02-sv-core/data/truth_svs_vaf{vaf}_{rep}.vcf",
    "ins":    "docs/benchmarking/02-sv-core/data/truth_ins_vaf{vaf}_{rep}.vcf",
    "invdel": "docs/benchmarking/02-sv-core/data/truth_invdel_vaf{vaf}_{rep}.vcf",
}

# How many configs to generate per type.
TYPE_COUNTS = {
    "snv":    900,
    "indel":  900,
    "sv":     400,
    "ins":    400,
    "invdel": 400,
}


def nearest_vaf_tag(vaf_fraction: float) -> int:
    """Return the closest available VAF tag (integer = vaf * 10000)."""
    vaf_int = vaf_fraction * 10000.0
    return min(AVAILABLE_VAF_TAGS, key=lambda t: abs(t - vaf_int))


def sample_params(rng: random.Random) -> dict:
    """Sample a random set of simulation parameters."""
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
    """Deterministic seed from config name, offset to avoid collisions with existing seeds."""
    return hash(name) % 100000 + 200000


def make_config_yaml(name: str, vtype: str, params: dict, truth_vcf: str, replicate: str) -> str:
    """Produce the varforge YAML config string."""
    vaf_tag = nearest_vaf_tag(params["vaf"])
    vaf_label = f"{vaf_tag / 10000:.4%}"  # e.g. "0.0500%"
    seed = config_seed(name)
    ref = REFERENCES[vtype]

    lines = [
        f"# varforge config: {name.upper()} — VAF {vaf_label}",
        f"reference: {ref}",
        "",
        "output:",
        f"  directory: bigdata/experiments/02-ml-single-strand/results/sim_{name}",
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


def make_params_json(name: str, vtype: str, params: dict, replicate: str) -> dict:
    """Produce the params.json dict alongside the config."""
    vaf_tag = nearest_vaf_tag(params["vaf"])
    return {
        "coverage": round(params["coverage"], 2),
        "family_size_mean": round(params["family_size_mean"], 3),
        "pcr_cycles": params["pcr_cycles"],
        "fragment_mean": round(params["fragment_mean"], 2),
        "variant_type": vtype,
        "vaf": round(params["vaf"], 6),
        "param_variant": "ml2",
        "replicate": replicate,
        "vaf_tag": vaf_tag,
    }


def generate_configs(seed: int = 42) -> list[str]:
    """Generate all 3000 configs. Returns a list of config names."""
    rng = random.Random(seed)
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
    names = []

    for vtype, count in TYPE_COUNTS.items():
        for idx in range(1, count + 1):
            # Sample params with a fresh draw.
            p = sample_params(rng)

            # Assign replicate randomly.
            replicate = rng.choice(["a", "b"])

            # VAF tag for the truth VCF filename.
            vaf_tag = nearest_vaf_tag(p["vaf"])
            vaf_tag_str = f"{vaf_tag:04d}"

            # Config name.
            vaf_label = int(round(p["vaf"] * 10000))
            name = f"{vtype}_ml2_vaf{vaf_label:04d}_{idx:04d}"

            # Truth VCF path.
            truth_vcf = TRUTH_VCF_TEMPLATES[vtype].format(
                vaf=vaf_tag_str, rep=replicate
            )

            # Write YAML.
            yaml_content = make_config_yaml(name, vtype, p, truth_vcf, replicate)
            yaml_path = CONFIGS_DIR / f"{name}.yaml"
            yaml_path.write_text(yaml_content)

            # Write params.json.
            params_data = make_params_json(name, vtype, p, replicate)
            params_path = CONFIGS_DIR / f"{name}_params.json"
            params_path.write_text(json.dumps(params_data, indent=2) + "\n")

            names.append(name)

    # Write batch manifest.
    manifest_path = CONFIGS_DIR / "ml2_manifest.json"
    manifest_path.write_text(json.dumps(names, indent=2) + "\n")

    return names


if __name__ == "__main__":
    names = generate_configs()
    print(f"Generated {len(names)} configs.")
    print(f"Manifest written to bigdata/experiments/02-ml-single-strand/configs/ml2_manifest.json")
