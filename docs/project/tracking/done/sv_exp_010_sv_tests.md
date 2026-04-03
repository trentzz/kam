# SV-EXP-010: Unit and Integration Tests for All New SV Types

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_002_invdel_classification.md, sv_exp_003_novel_insertion.md, sv_exp_006_fusion_classify.md, sv_exp_009_vcf_new_sv.md
**Status**: todo

## Goal

Write a comprehensive test suite covering all three new SV types end-to-end:
from synthetic path evidence through classification, threshold application, and
VCF output. Done looks like: a test file `kam-call/tests/sv_new_types.rs` (or
equivalent) with at minimum two tests per type (correct classification and
correct VCF output), all passing under `cargo test`.

## Steps

1. Read the existing SV test files in `kam-call/tests/` to understand the
   test patterns (how synthetic `PathEvidence` is constructed, how the caller
   is invoked).

2. For each new SV type, write:

   **InvDel**
   - `test_invdel_classification`: construct a path with an inverted deletion
     signature; assert `classify_variant` returns `InvDel`.
   - `test_invdel_vcf`: run the caller on this path; assert the VCF record
     contains `SVTYPE=INVDEL`.

   **NovelInsertion**
   - `test_novel_insertion_classification`: construct a path with a large
     insertion of non-repetitive sequence; assert `classify_variant` returns
     `NovelInsertion`.
   - `test_novel_insertion_vcf`: assert VCF record contains `SVTYPE=INS;NOVEL=1`.

   **Fusion**
   - `test_fusion_classification`: construct a path with k-mers annotated
     across two target IDs; assert `classify_variant` returns `Fusion`.
   - `test_fusion_vcf`: assert two BND VCF records are emitted with correct
     ALT notation.

3. Add regression tests:
   - A plain deletion is not reclassified as `InvDel`.
   - A tandem duplication is not reclassified as `NovelInsertion`.
   - A single-target path is not classified as `Fusion`.

4. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- Constructing synthetic `PathEvidence` may require helper builder functions.
  Add a `#[cfg(test)]` module with builders rather than duplicating
  construction code across tests.
- If the full pipeline integration test for fusions is complex, defer it to
  the BENCH-SVN-003 task. The unit tests here cover correctness at the
  classifier level.
