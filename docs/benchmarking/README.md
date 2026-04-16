# benchmarking

Benchmark definitions for kam. Each subdirectory covers one class of variants. Configs, scripts, reference data, and small result files are tracked here. Large generated files (varforge FASTQs, simulation outputs) live in `bigdata/benchmarking/` and are uploaded to Nextcloud.

## Subfolders

| Folder | Variant class | bigdata mirror |
|---|---|---|
| `01-snvindel/` | SNV and indel VAF sweep (0.05%–100%) | `bigdata/benchmarking/01-snvindel/` |
| `02-sv-core/` | Core SV types: large deletion, tandem duplication, inversion | `bigdata/benchmarking/02-sv-core/` |
| `03-sv-extended/` | Extended SV types: insertion, inversion-deletion, novel insertion | `bigdata/benchmarking/03-sv-extended/` |
| `04-comparison/` | Alignment-based baseline comparison | `bigdata/benchmarking/04-comparison/` |
| `05-public/` | Public datasets: COLO829, SEQC2, UMI benchmark | `bigdata/benchmarking/05-public/` |
| `06-runtime/` | Runtime profiling vs alignment-based pipeline | `bigdata/benchmarking/06-runtime/` |
| `07-snvindel-ml-boost-v1/` | ML model evaluation on real titration data | `bigdata/benchmarking/07-snvindel-ml-boost-v1/` |
| `resources/` | Reference data files shared across benchmarks | — |
| `scripts/` | Shared scoring and figure generation scripts | — |

## How to run

Each subfolder has its own README with dataset-specific instructions. The general pattern is:

1. Download or generate FASTQs into `bigdata/benchmarking/<folder>/`.
2. Run kam using the configs in that subfolder.
3. Score results using `scripts/score_all_results.py`.
4. Generate figures using `scripts/generate_paper_figures.py`.
