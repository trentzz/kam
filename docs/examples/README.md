# examples

Example configuration files for common kam use cases. Each file is a complete, runnable TOML config. Copy the file that matches your chemistry and workflow, then adjust paths and parameters for your sample.

## Contents

| File | Chemistry | Workflow |
|---|---|---|
| `twist-umi-duplex.toml` | Twist UMI duplex (5 bp inline UMI) | Full duplex variant calling |
| `tumour-informed.toml` | Twist UMI duplex | Tumour-informed mode with known variant targets |
| `ctdna-monitoring.toml` | Twist UMI duplex | ctDNA monitoring against a fixed variant list |
| `discovery-mode.toml` | Twist UMI duplex | De novo discovery mode |
| `fast-mode.toml` | Twist UMI duplex | Fast mode for screening |
| `high-sensitivity.toml` | Twist UMI duplex | High-sensitivity settings for very low VAF |
| `high-specificity.toml` | Twist UMI duplex | Conservative settings to minimise false positives |
| `large-panel.toml` | Twist UMI duplex | Settings for panels >500 targets |
| `low-input.toml` | Twist UMI duplex | Low DNA input (low family size) |
| `minimal.toml` | Any | Minimal config with only required fields |
| `no-umi.toml` | No UMI | Standard paired-end without UMI |
| `research-discovery.toml` | Twist UMI duplex | Broad discovery for research use |
| `simplex-umi-8bp.toml` | Simplex 8 bp UMI | |
| `simplex-umi-9bp.toml` | Simplex 9 bp UMI | |
| `simplex-umi-12bp.toml` | Simplex 12 bp UMI | |
| `sv-detection.toml` | Twist UMI duplex | Structural variant detection |
| `sv-high-sensitivity.toml` | Twist UMI duplex | High-sensitivity SV detection |
| `thorough-mode.toml` | Twist UMI duplex | Thorough mode: slower but more complete |
| `fusion-detection.toml` | Twist UMI duplex | Fusion gene detection |
| `ti-pipeline/` | Twist UMI duplex | Tumour-informed pipeline configs |
