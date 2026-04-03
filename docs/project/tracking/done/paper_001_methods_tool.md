# PAPER-001: Methods Section — Tool Description

**Epic**: PAPER (docs/claudetracking/overallplans/PAPER.md)
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Write the tool description subsection of the methods: kam's architecture,
the config system, and multi-chemistry support. Done looks like: a complete
draft (≈600 words) committed to `docs/paper/sections/methods_tool.md`, with
accurate technical detail matching the implemented code.

## Steps

1. Read `docs/planning/rust_workspace_architecture.md` and
   `docs/planning/core_data_model.md` to ensure the description is accurate.

2. Write `docs/paper/sections/methods_tool.md` covering:

   **Architecture (≈200 words)**
   - Five-stage pipeline: assemble, index, pathfind, call, output.
   - Each stage as a separate Rust crate.
   - Molecule as the atomic unit of evidence throughout.
   - No alignment step: k-mer-based detection.

   **Configuration (≈150 words)**
   - TOML config file system (CONFIG-TOML epic).
   - Layered resolution: defaults → config file → CLI overrides.
   - Example: `kam run --config twist-umi-duplex.toml --r1 R1.fq.gz`.

   **Chemistry support (≈150 words)**
   - Generalised UMI length and skip length (CHEM-CONFIG epic).
   - Named presets for Twist duplex, simplex 12 bp, simplex 9 bp.
   - Duplex UMI canonical pairing: `min(ZA+ZB, ZB+ZA)`.

   **Variant calling (≈100 words)**
   - Binomial confidence model.
   - SV-specific evidence model (mean_variant_specific_molecules).
   - Tumour-informed mode.

3. Write in active voice, Australian English, no em dashes, no padding.

4. Note any sections that require results to be filled in (e.g. performance
   numbers) with `[TODO: fill after benchmarks]` placeholders.

## Notes

- Do not describe the statistical model in detail here; that goes in the
  benchmarking methods subsection (paper_002).
- Target audience: computational biologists familiar with bioinformatics
  tools but not necessarily with Rust or k-mer methods.
- This section can be drafted before benchmarks are complete.
