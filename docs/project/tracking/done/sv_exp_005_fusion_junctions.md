# SV-EXP-005: Fusion Junction K-mer Generation

**Epic**: SV-EXPAND (docs/claudetracking/overallplans/SV-EXPAND.md)
**Priority**: high
**Depends on**: sv_exp_004_fusion_design.md
**Status**: todo

## Goal

Implement a utility (config file or CLI sub-command) that generates junction
k-mers for specified fusion partner pairs. Done looks like: given a partner-pair
specification (partner A sequence + partner B sequence), the utility writes a
FASTA of `k - 1` junction k-mers that span the A-B boundary and can be used
as targets for the k-mer index.

## Steps

1. Read `docs/research/fusion_detection_design.md` (written in SV-EXP-004).
   Follow the design decisions documented there.

2. Define the partner-pair input format as decided in SV-EXP-004. Likely a
   FASTA file with two sequences per fusion event, or a TSV with columns
   `name`, `partner_a_seq`, `partner_b_seq`.

3. Implement `generate_fusion_junctions(partner_a: &[u8], partner_b: &[u8],
   k: usize) -> Vec<Vec<u8>>` in a new file `kam-pathfind/src/fusion.rs` (or
   in `kam/src/fusion_prep.rs` if it is a pre-processing step):
   - Concatenate last `k-1` bases of A with first `k-1` bases of B.
   - Extract all k-mers from the junction string (there are `k - 1` of them).
   - Return as a `Vec<Vec<u8>>`.

4. Annotate each junction k-mer with both `target_id` values (A and B) so the
   classifier can later identify them as spanning two targets.

5. Write unit tests:
   - Correct k-mer count for given k.
   - All k-mers contain exactly the A-B boundary.
   - K-mers that appear in the normal reference are flagged (or at minimum:
     a test that detects the overlap).

6. Wire into the pipeline: junction k-mers are added to the allowlist (or a
   separate junction index) before the main `kam index` step.

7. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- The exact wiring into the index depends on decisions in SV-EXP-004.
  Do not make architectural decisions here; implement what the design doc says.
- If the design doc deferred the junction k-mer approach, do not implement it
  here. Instead, update the task to match the agreed design.
