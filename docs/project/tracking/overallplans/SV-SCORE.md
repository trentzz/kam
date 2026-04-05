# SV-SCORE: SV Junction Evidence Model

**Status**: complete
**Priority**: high
**Branch**: epic/SV-SCORE (merged to main 2026-03-24)

## Goal

Fix the `min_molecules` bottleneck that makes large SV paths appear to have 1–3
supporting molecules regardless of true VAF.

## Background

Large SVs (inversions, insertions, tandem duplications ≥50bp) are detected by
walking a path through the de Bruijn graph. The path spans the entire SV region:
~70–150 k-mers for a 100bp event. The caller uses `min_molecules` across all
k-mers to estimate the alt-allele evidence count.

For an inversion path, each k-mer inside the inverted region is covered only by
molecules that span that exact position within the inversion. With 5000x coverage
at 1% VAF (~50 mutant molecules), the worst-covered k-mer deep in the inversion
has 1–3 supporting molecules. This sets `min_molecules = 1–3` regardless of the
true VAF, making the reported NALT always low and the VAF estimate always far below
truth.

Evidence (from `docs/research/sv_large_detection_investigation.md`):
- At 10% VAF, NALT=41 for the INV while the DUP (short path) has NALT ≫ 100
- Reported INV VAF is always 10–40% of truth VAF due to the min bottleneck
- The DUP (1bp junction path) has no bottleneck and reports inflated VAF

## Design

### Core change: mean_variant_specific_molecules

Add a new field to `PathEvidence` in `kam-pathfind/src/score.rs`:

```rust
/// Mean `n_molecules` across variant-specific k-mers only.
///
/// For large SV paths, this is a better evidence estimate than min_molecules
/// because the path traverses ~100+ k-mers, and the minimum over that many
/// positions bottlenecks at 1–3 even when true VAF is high.
///
/// Variant-specific k-mers are those not shared with the reference path.
/// Set to 0.0 for the reference path. Computed by score_and_rank_paths.
pub mean_variant_specific_molecules: f32,
```

Computed in `score_and_rank_paths` using the existing `ref_kmer_set` logic (same
approach as `min_variant_specific_duplex`).

### Caller change: SV-specific evidence

In `kam-call/src/caller.rs`, `call_variant`, use `mean_variant_specific_molecules`
for `k` and `total` when the variant type is SV:

```rust
let k = if is_sv_type(variant_type) {
    alt_evidence.mean_variant_specific_molecules.round() as u32
} else {
    alt_evidence.min_molecules
};
```

This gives a representative molecule count for large SV paths without inflating
FPs from error k-mers (which have min_molecules=1 and mean_variant_specific=0.0).

### CLI change: no new flags

This change is internal to the scoring model. No new CLI flags required.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-SCORE-001 | done/sv_score_001_mean_vs_molecules.md | done |
| SV-SCORE-002 | done/sv_score_002_caller_sv_evidence.md | done |
| SV-SCORE-003 | done/sv_score_003_benchmark_rerun.md | done |

## Scope

- `kam-pathfind/src/score.rs` — add `mean_variant_specific_molecules` field and computation
- `kam-call/src/caller.rs` — use new field for SV types
- `kam-pathfind/src/` and `kam-call/src/` tests — update for new field
- `docs/benchmarking/sv/` — re-run and re-score SV suite

## Out of scope

- Changes to `kam-core` types (requires explicit approval)
- Changes to the de Bruijn graph or path walking
- Strand bias behaviour (tracked in SV-STRAND epic)
