# PERF-001: Replace Vec::contains with HashSet in assembler clustering

**Epic**: none (standalone fix)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Replace `Vec::contains` (O(n) per lookup) with `HashSet::contains` (O(1)) in
the UMI family clustering loop in `kam-assemble/src/assembler.rs`. Done looks
like: the clustering loop uses a `HashSet`, and a benchmark or comment
documents the improvement.

## Steps

1. Read `kam-assemble/src/assembler.rs`, focusing on `cluster_read_indices`
   or the equivalent clustering function (around line 492 for family grouping).

2. Identify all `Vec::contains` calls inside the clustering loop. Replace the
   backing container with a `HashSet` where appropriate.

3. Also fix finding OLD-1.3 in the same pass: the family size cast to `u8`:
   ```rust
   let n_fwd = reads.iter().filter(|r| r.is_forward).count();
   ```
   If `n_fwd > 255`, this silently wraps. Use `u8::try_from(n_fwd).unwrap_or(u8::MAX)`
   or saturating cast.

4. Move this task file to done and update Status.

5. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review findings OLD-3.1 (Vec::contains) and OLD-1.3 (u8 cast).
- Both fixes are in the same file and are small enough to combine.
