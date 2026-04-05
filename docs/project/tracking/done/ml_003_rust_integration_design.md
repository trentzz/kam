# ML-003: Design Rust integration for ML model inference

**Epic**: ML-BOOST (overallplans/ML-BOOST.md)
**Priority**: medium
**Depends on**: ML-002
**Status**: todo

## Goal

Design the concrete integration plan for running LightGBM inference inside
the kam Rust binary. Write a design document (not code) covering CLI flag
design, model loading, inference at call time, and output format changes.

## Steps

1. Research available options:
   - `lightgbm` Rust crate (wraps the C API)
   - `ort` (ONNX Runtime Rust) — LightGBM can export to ONNX
   - Python subprocess at call time (simplest, no Rust changes to inference)
2. Evaluate each option: build complexity, binary size, performance, portability.
3. Design the CLI flag: `--ml-model <path>` on `kam call`. Optional. If not
   provided, behaviour is identical to current.
4. Design the output change: add `ml_prob` (f32) and `ml_filter` (string:
   PASS or MLFiltered) columns to the call TSV.
5. Write the design document to
   `docs/planning/ml_integration_design.md`.
6. Include a concrete build plan: which crate to add, which structs change,
   which functions to add.

## Notes

- Do not implement any Rust code in this task. Design only.
- The design must not require changes to kam-core types without explicit
  approval (per CLAUDE.md).
- Consider the `lightgbm-rs` crate as a first option.

## Success criteria

- [ ] `docs/planning/ml_integration_design.md` exists
- [ ] Document covers CLI design, model loading, inference, output format
- [ ] At least two Rust integration options compared with tradeoffs
- [ ] Concrete recommended approach stated with rationale
