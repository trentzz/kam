# Experiments

## 1. Titration Comparison

**Goal:** Validate tumour-informed sensitivity versus the alignment-based baseline on all titration samples.

- Input: all 24 titration samples (5/15/30 ng input × 8 VAF levels).
- Mode: tumour-informed.
- Baseline: RaSCALL results from `/mnt/tzeng-local/tzeng-thesis/`.
- Analysis: per-variant concordance. For every variant, record whether kam and the alignment-based baseline detect it, at what VAF, and with what confidence score.
- Output: sensitivity table, precision table, per-variant concordance table.

## 2. SV Expansion Benchmark

**Goal:** Demonstrate SV detection across all supported types.

- SV types: fusions, translocations, large deletions, inversions, duplications.
- VAF range: 0.5% to 5%.
- Input: varforge simulations.
- Output: sensitivity heatmap by SV type and VAF level.

## 3. Public Dataset Benchmarks

**Goal:** Show that kam generalises beyond the titration set.

- Select diverse public datasets with known ground truth.
- Run kam in tumour-informed mode.
- Compare against published alignment-based results.
- Output: sensitivity and precision on each dataset.

## 4. Cross-Chemistry Test

**Goal:** Demonstrate that kam is not specific to Twist chemistry.

- Run kam with different `config.toml` chemistry settings on non-Twist UMI data.
- Report sensitivity and precision on non-Twist samples.
- Output: per-chemistry results table.

## 5. Runtime Benchmark

**Goal:** Quantify the speed advantage over alignment-based pipelines.

- Run kam and the alignment-based pipeline on matched input datasets.
- Measure wall-clock time for each stage.
- Output: runtime comparison bar chart.

## 6. Per-Variant Deep Dive

**Goal:** Build a complete picture of where each method succeeds and fails.

- For every variant in the titration set, document:
  - Detection status in kam (yes/no, VAF threshold, confidence).
  - Detection status in alignment-based (yes/no, VAF threshold, confidence).
  - Discordance analysis: why does one method succeed where the other fails?
- Output: annotated per-variant table for supplementary material, summary in main text.
