# kam — Alignment-Free Variant Detection for Duplex UMI Sequencing

## Project Description

kam is a fully alignment-free variant detection pipeline for duplex UMI sequencing (Twist Biosciences chemistry), built in Rust. It replaces the current pipeline of HUMID → Jellyfish → km → kmtools with a unified Rust workspace that preserves molecule-level information throughout.

Target use case: detecting somatic variants at 0.1% VAF in liquid biopsy ctDNA panel samples.

## Chemistry: Twist UMI Duplex

- Read structure: `5M2S+T` on both R1 and R2
- 5bp random UMI (1,024 possible values — collision-aware design required)
- 2bp skip/spacer (monotemplate, QC signal)
- Duplex identification: canonical pair = `min(ZA+ZB, ZB+ZA)` lexicographically
- See `docs/research/twist_umi_chemistry.md` for full details

## Architecture Decisions — DO NOT CHANGE WITHOUT ASKING

- **Molecule is the atomic unit**, not Read. All counts, grouping, and evidence are molecule-level.
- **kam-core defines shared types** — all other crates depend on it. Do not modify kam-core types without explicit approval (commit message must contain "APPROVED").
- **KmerIndex stores MoleculeEvidence**, not integer counts. This is the fundamental advantage over Jellyfish.
- **All RNG must use seeded deterministic generators** for Nextflow cache compatibility.
- **No unsafe code** without explicit comment justification.
- **All public functions must have doc comments** with examples.
- **Targeted mode first**, de novo discovery second.

## Workspace Crates

```
kam-core     — shared types (Molecule, ConsensusRead, MoleculeEvidence, traits)
kam-assemble — molecule assembly from raw FASTQ (replaces HUMID)
kam-index    — k-mer indexing with molecule provenance (replaces Jellyfish)
kam-pathfind — de Bruijn graph construction and path walking (replaces km)
kam-call     — statistical variant calling with confidence intervals
kam          — integrated binary CLI
```

See `docs/project/devmanual/rust_workspace_architecture.md` for full architecture.

## Code Quality Rules

- `cargo test` must pass before marking any task done
- `cargo clippy -- -D warnings` must pass
- No `unwrap()` in library code — use proper error handling with `thiserror`
- Every function has a unit test
- Every module has an integration test using synthetic data
- Deterministic output for given inputs (no HashMap iteration order in output)

## Key Dependencies

- `needletail` — streaming FASTQ parsing
- `dashmap` — concurrent HashMap
- `rayon` — parallelism
- `serde` + `bincode` — serialisation
- `clap` — CLI
- `statrs` — statistical distributions
- `thiserror` — error types

## Documentation Pointers

### Research (prior art and background)
- `docs/research/humid_analysis.md` — HUMID algorithm, capabilities, and gaps
- `docs/research/tool_landscape.md` — UMI-tools, fgbio, UMICollapse, Sentieon, Jellyfish, km
- `docs/research/gaps_to_fill.md` — 7 differentiated value areas for kam
- `docs/research/twist_umi_chemistry.md` — Twist UMI chemistry details, parsing, grouping
- `docs/research/kmer_memory_strategies.md` — allowlist, Bloom filter, compressed representations
- `docs/research/graph_building_challenges.md` — anchor contamination, bridge problem, de novo vs targeted
- `docs/research/endpoint_fingerprinting.md` — UMI collision detection via template endpoints, bases/threshold tradeoffs
- `docs/research/twist_skip_bases.md` — skip base identity research, auto-detection approach, QC usage
- `docs/research/consensus_calling_algorithms.md` — majority vote, quality-weighted, Bayesian, duplex crossing, edge cases
- `docs/research/statistical_calling_models.md` — binomial/beta-binomial models, duplex-aware scoring, strand bias, background error model, validation strategy
- `docs/research/streaming_molecule_assembly.md` — sort-then-group, hash-partition, memory bounds, interaction with Hamming clustering

### Developer manual (architecture and design)
- `docs/project/devmanual/rust_workspace_architecture.md` — crate layout, dependencies, build order
- `docs/project/devmanual/core_data_model.md` — Molecule, ConsensusRead, MoleculeEvidence, traits
- `docs/project/devmanual/design_principles.md` — transparency, granular control, extreme logging, explain-why, no black boxes, reproducibility
- `docs/project/devmanual/logging_architecture.md` — per-file togglable logs, drop log format, zero-cost-when-off pattern, Nextflow interaction
- `docs/project/devmanual/nextflow_integration.md` — HPC config, QC JSON, morning report
- `docs/project/devmanual/development_workflow.md` — Claude Code loop, task format, safety boundaries
- `docs/project/devmanual/current_pipeline_analysis.md` — information loss analysis at each pipeline stage
- `docs/project/devmanual/de_novo_discovery_design.md` — 4 approaches explored, panel-aware de novo as near-term extension, deferred to Phase 2
- `docs/project/devmanual/output_format_specs.md` — annotated FASTQ tag definitions, QC JSON schemas per stage, TSV/VCF/JSON variant output, collision probability reporting
- `docs/project/devmanual/nextcloud.md` — Nextcloud download and upload instructions, bigdata/ structure

### Other
- `docs/project/features/` — feature specs (todo/inprogress/done)
- `docs/project/tracking/` — kanban task board
- `docs/project/investigations/` — diagnostic investigation writeups
- `docs/project/experiments/` — research experiments (ML, etc.)
- `docs/benchmarking/` — benchmark definitions, configs, scripts
- `docs/paper/` — publication materials
- `docs/manual/` — end-user documentation
- `docs/examples/` — example configuration files

## Guides

Detailed instructions for working in this repo live in `docs/project/devmanual/`.
Read the relevant guide before starting work in that area.

| Guide | When to read |
|-------|-------------|
| `docs/project/devmanual/benchmarking.md` | Running benchmarks, adding datasets, scoring results |
| `docs/project/devmanual/task-tracking.md` | Picking up tasks, creating epics, committing work |
| `docs/project/devmanual/investigation-docs.md` | Writing up a diagnostic investigation |
| `docs/claudeguide/release.md` | Cutting a release, building the binary, publishing to GitHub Releases |

## Autonomous Session Behaviour

### What To Do When Stuck
1. Design decision needed → write `// DESIGN_QUESTION: [question]` comment and move to next task
2. Tests fail after 3 attempts → write `// TODO: [what you tried]` and move on
3. **Never change interface definitions in kam-core** without explicit instruction
4. If unsure about scope, err on the side of smaller

### Current Sprint Goal
Benchmarking: varforge VAF sweeps for SNV/indel (0.5–5%) and SV new types
(INS, large DEL, INVDEL at 0.5–5%).
