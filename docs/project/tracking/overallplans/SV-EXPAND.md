# SV-EXPAND: New SV Types (Fusion, InvDel, NovelInsertion)

**Status**: done
**Priority**: high
**Branch**: epic/SV-EXPAND

## Goal

Extend kam's variant type system with three new structural variant classes:
`Fusion`, `InvDel`, and `NovelInsertion`. Update classification logic, path
walking allele extraction, threshold application, and VCF output. After this
epic, kam can classify and report all five SV types relevant to ctDNA liquid
biopsy panels.

## Motivation

The current `VariantType` enum covers `LargeDeletion`, `TandemDuplication`,
and `Inversion`. Clinically relevant SVs also include:

- **InvDel**: an inverted sequence with flanking deletions, common in
  chromosome rearrangements.
- **NovelInsertion**: insertion of sequence with no tandem homology to the
  reference, e.g. retrotransposon or viral integration.
- **Fusion**: junction spanning two distinct genomic loci (translocation or
  gene fusion), the most clinically actionable SV in haematological malignancies.

## Design

### Fusion detection

Fusion requires a separate k-mer generation step: build junction k-mers that
span the boundary of two target regions (partner A last `k-1` bp + partner B
first `k-1` bp). A path is classified as a Fusion when its k-mers resolve to
two distinct `target_id` values. This requires extending the target FASTA
concept and the classification logic. A design document (sv_exp_004) is
required before implementation of sv_exp_005 and sv_exp_006.

### InvDel classification

Re-classify `LargeDeletion` calls as `InvDel` when the differing region
between ref and alt paths contains a reverse-complement segment. The classifier
already has the ref and alt sequences; add an RC subsequence check.

### NovelInsertion classification

Re-classify `TandemDuplication` calls as `NovelInsertion` when the inserted
sequence (≥50 bp) does not repeat a nearby reference window. Short tandem
repeats and homopolymers are excluded from this re-classification.

### Thresholds

Apply SV-specific thresholds (`sv_min_confidence`, `sv_min_alt_molecules`,
`sv_strand_bias_threshold`) to all three new types using the existing
`is_sv_type` predicate.

### VCF output

- `InvDel`: `SVTYPE=INVDEL`
- `NovelInsertion`: `SVTYPE=INS` with `NOVEL=1` INFO flag
- `Fusion`: BND notation (two VCF records per event, standard PBSV/SVABA
  convention)

## Child tasks

| ID | File | Status |
|----|------|--------|
| SV-EXP-001 | done/sv_exp_001_new_variant_types.md | done |
| SV-EXP-002 | done/sv_exp_002_invdel_classification.md | done |
| SV-EXP-003 | done/sv_exp_003_novel_insertion.md | done |
| SV-EXP-004 | done/sv_exp_004_fusion_design.md | done |
| SV-EXP-005 | done/sv_exp_005_fusion_junctions.md | done |
| SV-EXP-006 | done/sv_exp_006_fusion_classify.md | done |
| SV-EXP-007 | done/sv_exp_007_multi_coord_allele.md | done |
| SV-EXP-008 | done/sv_exp_008_sv_thresholds.md | done |
| SV-EXP-009 | done/sv_exp_009_vcf_new_sv.md | done |
| SV-EXP-010 | done/sv_exp_010_sv_tests.md | done |

## Scope

- `kam-core/src/variant.rs` — add new `VariantType` variants
- `kam-call/src/classify.rs` — InvDel and NovelInsertion classification
- `kam-call/src/caller.rs` — threshold routing for new types
- `kam-call/src/output.rs` or `kam/src/output.rs` — VCF BND and SVTYPE output
- `kam-pathfind/src/` — fusion junction k-mer generation utility
- `docs/research/` — fusion design document

## Out of scope

- Changes to de Bruijn graph construction
- Repeat annotation databases (use local sequence context only)
- Somatic vs germline classification
