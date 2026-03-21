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

See `docs/planning/rust_workspace_architecture.md` for full architecture.

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

### Research (background context)
- `docs/research/humid_analysis.md` — HUMID algorithm, capabilities, and gaps
- `docs/research/tool_landscape.md` — UMI-tools, fgbio, UMICollapse, Sentieon, Jellyfish, km
- `docs/research/gaps_to_fill.md` — 7 differentiated value areas for kam
- `docs/research/twist_umi_chemistry.md` — Twist UMI chemistry details, parsing, grouping
- `docs/research/kmer_memory_strategies.md` — allowlist, Bloom filter, compressed representations
- `docs/research/graph_building_challenges.md` — anchor contamination, bridge problem, de novo vs targeted
- `docs/research/endpoint_fingerprinting.md` — UMI collision detection via template endpoints, bases/threshold tradeoffs
- `docs/research/twist_skip_bases.md` — skip base identity research, auto-detection approach, QC usage
- `docs/research/consensus_calling_algorithms.md` — majority vote, quality-weighted, Bayesian, duplex crossing, edge cases
- `docs/research/de_novo_discovery_design.md` — 4 approaches explored, panel-aware de novo as near-term extension, deferred to Phase 2
- `docs/research/statistical_calling_models.md` — binomial/beta-binomial models, duplex-aware scoring, strand bias, background error model, validation strategy
- `docs/research/streaming_molecule_assembly.md` — sort-then-group, hash-partition, memory bounds, interaction with Hamming clustering
- `docs/research/output_format_specs.md` — annotated FASTQ tag definitions, QC JSON schemas per stage, TSV/VCF/JSON variant output, collision probability reporting

### Planning (architecture and design)
- `docs/planning/current_pipeline_analysis.md` — information loss analysis at each pipeline stage
- `docs/planning/rust_workspace_architecture.md` — crate layout, dependencies, build order
- `docs/planning/core_data_model.md` — Molecule, ConsensusRead, MoleculeEvidence, traits
- `docs/planning/nextflow_integration.md` — HPC config, QC JSON, morning report
- `docs/planning/design_principles.md` — transparency, granular control, extreme logging, explain-why, no black boxes, reproducibility
- `docs/planning/logging_architecture.md` — per-file togglable logs, drop log format, zero-cost-when-off pattern, Nextflow interaction
- `docs/planning/development_workflow.md` — Claude Code loop, task format, safety boundaries

### Other
- `docs/features/` — feature specs (todo/inprogress/done)
- `docs/benchmarking/` — performance benchmarks, sensitivity metrics
- `docs/paper/` — publication materials
- `docs/manual/` — end-user documentation

## Investigation Documentation

Whenever a diagnostic investigation uncovers a non-obvious finding (root cause of a bug,
unexpected performance characteristic, surprising benchmark result, failed hypothesis), write
it up in `docs/research/` immediately. Use the pattern established in
`sensitivity_investigation.md` and `anchor_missing_investigation.md`:

1. State the symptom and the diagnostic data that revealed it.
2. State the hypothesis tested and why it was right or wrong.
3. State the actual root cause with supporting evidence.
4. State what would actually fix it and what the expected improvement is.
5. State what was implemented and what the measured result was.

This log is essential context for future sessions and paper writing. Do not just fix the code
and move on — document the reasoning.

---

## Autonomous Session Behaviour

### What To Do When Stuck
1. Design decision needed → write `// DESIGN_QUESTION: [question]` comment and move to next task
2. Tests fail after 3 attempts → write `// TODO: [what you tried]` and move on
3. **Never change interface definitions in kam-core** without explicit instruction
4. If unsure about scope, err on the side of smaller

### Task Queue
- `docs/claudeloop/queue/` — tasks to pick up (sorted by filename)
- `docs/claudeloop/in_progress/` — currently being worked on
- `docs/claudeloop/done/` — completed tasks
- `docs/claudeloop/needs_review/` — tasks Claude got stuck on
- `docs/claudeloop/open_questions.md` — design decisions needing user input

### Current Sprint Goal
Implement kam-core types and kam-assemble parser for Twist UMI duplex data.
