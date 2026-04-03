# Open Design Questions

Questions needing user input before implementation can proceed.

## Q1: Minimum Template Length — RESOLVED

Parser is k-agnostic. No hardcoded minimum. User can optionally set `--min-template-length` but default is to parse everything and let downstream stages skip reads too short for the chosen k. Dropped reads are logged with reason when logging enabled.

See `docs/planning/logging_architecture.md` for the broader logging design.

## Q2: UMI Quality Threshold Default — RESOLVED

Default Q20 per individual UMI base (any of the 10 UMI bases across R1+R2 below Q20 → read dropped). Configurable via `--min-umi-quality`. At Q20, ~10% of reads have at least one UMI error across 10 bases — Hamming distance clustering (distance 1) handles most remaining errors. Dropped reads logged with which base failed and its quality value.

## Q3: Endpoint Fingerprint Parameters — RESOLVED

Default: 8 bases per endpoint, 4 bit difference threshold. Configurable via `--fingerprint-bases` and `--fingerprint-max-diff`. 8bp exactly fills a u64 (no hashing), 4-bit-diff tolerates ~2 sequencing errors across endpoints while catching genuine collisions. See `docs/research/endpoint_fingerprinting.md` for full tradeoff analysis.

## Q4: Expected Skip Bases for Twist Chemistry — RESOLVED

The actual 2bp spacer sequence is **not publicly documented** by Twist. It's a fixed monotemplate consistent within an adapter lot, but cross-lot consistency is unconfirmed. No existing tools (fgbio, fastp, etc.) validate or hardcode it.

Approach: **auto-detect from data** — parse positions 5-6 across reads, the dominant dinucleotide (>95% frequency) is the skip sequence for that lot. Report distribution in QC. Optionally validate against `--expected-skip-bases <seq>` if user knows their lot's sequence. Flag mismatches in QC and logs but don't hard-reject by default. See `docs/research/twist_skip_bases.md`.

To get a definitive answer on lot-to-lot consistency, contact Twist technical support or check the adapter Certificate of Analysis.

## Q5: Binary Serialization Format for Inter-Stage Data — RESOLVED

`bincode` for the hot path between pipeline stages (maximum performance). Users can export to TSV, CSV, or JSON via `kam export` for inspection or cross-language access. Annotated FASTQ is the molecule handoff format between stages (human-readable, tool-compatible); bincode is used for the k-mer index and internal checkpoints.

Additionally: `kam explore` will be an interactive REPL-style CLI that reads bincode files directly for searching, filtering, and selective export without requiring full dump. See `docs/features/todo/interactive_explorer.md`.

Bincode files should include a header with file type, version, record count, and summary stats so the explorer can show `summary` instantly.

## Q6: Consensus Calling Algorithm — RESOLVED

Phase 1: quality-weighted majority vote (fast, simple, testable). Phase 2: Bayesian posterior with log-space computation (calibrated confidence scores). Duplex crossing: pick higher-confidence SSC at disagreement positions, flag as non-duplex. `--strict-duplex` flag masks disagreements as N instead. All parameters configurable. See `docs/research/consensus_calling_algorithms.md` for full analysis.

## Q7: De Novo Discovery Priority — RESOLVED

Targeted detection first. De novo deferred to later phase. Targeted mode supports:
- Specific SNVs/indels (exact variant to check for)
- Regions (find any variants within a defined region — a natural near-term extension)

See `docs/research/de_novo_discovery_design.md` for full exploration of 4 de novo approaches, with panel-aware de novo (Approach 1) as the recommended near-term extension and two-pass molecule-aware counting (Approach 4) as the medium-term path.

## Q8: Output Formats — RESOLVED

Support all formats, let the user pick what fits their workflow: TSV (default, km-compatible), CSV, JSON, VCF. Multiple can be emitted simultaneously via `--output-format tsv,json`. Same data, different serialization — minimal implementation cost. See `docs/research/output_format_specs.md` for column/field definitions.

## Q9: Multi-k Support Priority — RESOLVED

Single k first (default 31), multi-k as a later enhancement. User wants this as a feature eventually — variants called at multiple k values with agreement required would be higher confidence. The `KmerIndex` trait and `VariantGraph` trait already support instantiating multiple instances at different k values, so the architecture doesn't need to change — it's an integration task on top of the existing design.
