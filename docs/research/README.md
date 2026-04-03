# research

Background reading and prior art. These documents capture the analysis done before and during design of kam. They are reference material, not instructions.

## Contents

| File | Topic |
|---|---|
| `humid_analysis.md` | HUMID algorithm: capabilities, limitations, and gaps |
| `tool_landscape.md` | UMI-tools, fgbio, UMICollapse, Sentieon, Jellyfish, km |
| `gaps_to_fill.md` | 7 differentiated value areas where kam improves on existing tools |
| `twist_umi_chemistry.md` | Twist UMI duplex chemistry: read structure, parsing, grouping |
| `kmer_memory_strategies.md` | Allowlist, Bloom filter, and compressed k-mer representations |
| `graph_building_challenges.md` | Anchor contamination, bridge problem, de novo vs targeted |
| `endpoint_fingerprinting.md` | UMI collision detection via template endpoints |
| `twist_skip_bases.md` | Skip base identity research and auto-detection approach |
| `consensus_calling_algorithms.md` | Majority vote, quality-weighted, Bayesian, and duplex crossing |
| `statistical_calling_models.md` | Binomial and beta-binomial models, duplex-aware scoring, strand bias |
| `streaming_molecule_assembly.md` | Sort-then-group, hash-partition, memory bounds |

Investigation documents (diagnoses of specific problems encountered during development) live in `project/investigations/`. Internal design documents (implementation decisions) live in `project/devmanual/`.
