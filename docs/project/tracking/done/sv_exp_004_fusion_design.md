# SV-EXP-004: Fusion Detection Design Document

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_001_new_variant_types.md
**Status**: todo

## Goal

Write a design document for fusion/translocation detection in kam. The document
must answer the key architectural questions before any code is written for
SV-EXP-005 and SV-EXP-006. Done looks like: a committed Markdown file at
`docs/research/fusion_detection_design.md` that covers junction k-mer
generation, graph representation, classification, and known failure modes.

## Steps

1. Read:
   - `docs/research/graph_building_challenges.md`
   - `docs/research/endpoint_fingerprinting.md`
   - `docs/planning/core_data_model.md`
   - `kam-pathfind/src/` (graph construction and path walking) for current
     single-target paradigm.

2. Answer each of these questions in the document:

   **Junction k-mer generation**
   - A fusion between target A (e.g. BCR exon 14) and target B (e.g. ABL1
     exon 2) produces a junction sequence: `...last_k-1_bases_of_A +
     first_k-1_bases_of_B...`.
   - How are fusion partner pairs specified by the user? FASTA of breakpoint
     sequences? A partner-pair TSV?
   - How many junction k-mers per pair? (Answer: `k - 1` overlapping k-mers.)

   **Graph representation**
   - Currently, each target has its own de Bruijn graph. For fusions, do we
     build a joint graph across both targets, or a separate junction graph?
   - What happens when junction k-mers appear in the graph of only one target?

   **Classification**
   - How does the classifier know a path spans two targets? Via `target_id`
     annotations on k-mers?
   - What is the minimum evidence (molecules) for a fusion call?

   **Failure modes**
   - Paralogous sequences (both partners appear in the same genomic window).
   - Repeated junction k-mers (k-mer also occurs in normal sequence).
   - Unknown partner (de novo fusion discovery, out of scope).

3. Include a proposed implementation plan: which crates change, in what order,
   and what the interfaces look like.

4. Write the document to `docs/research/fusion_detection_design.md`.

## Notes

- This is a design task, not a code task. No Rust changes.
- The document informs SV-EXP-005 (junction FASTA generation) and SV-EXP-006
  (classification logic). Both tasks should be blocked until this document is
  reviewed and approved.
- If the design reveals that fusion detection requires significant architectural
  changes beyond the current sprint, document that clearly and flag for
  re-prioritisation.
