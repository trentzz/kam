# SV-CLASSIFY: Fix NovelInsertion Misclassification

**Status**: todo
**Priority**: high
**Branch**: epic/SV-CLASSIFY

## Goal

Diagnose and fix the bug where NovelInsertion is classified as TandemDuplication
by the `is_tandem_duplication()` function. Verify the fix with unit tests and the
full benchmark suite.

## Motivation

Current benchmark results show all NovelInsertion variants reported as
TandemDuplication (SVTYPE=DUP instead of SVTYPE=INS). This undermines the
classification system added in SV-EXPAND and makes NovelInsertion detection
results meaningless.

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-CLS-001 | todo/sv_cls_001_investigate.md | todo |
| SV-CLS-002 | todo/sv_cls_002_fix.md | todo |
| SV-CLS-003 | todo/sv_cls_003_regression_tests.md | todo |
| SV-CLS-004 | todo/sv_cls_004_rerun_novins.md | todo |
| SV-CLS-005 | todo/sv_cls_005_compare.md | todo |

## Scope

- `kam-call/src/caller.rs` — classification logic fix
- `docs/benchmarking/sv_new/results/kam_novins_*` — re-run results
- `docs/research/` — investigation document
