# BUG-003: call_ssc_bases uses template_r1 for both strands

**Epic**: none (standalone bug fix)
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

Fix `call_ssc_bases` in `kam-assemble/src/assembler.rs` to use `template_r2`
for reverse-strand reads instead of always using `template_r1`. Done looks
like: the duplex consensus uses R1 templates for the forward SSC and R2
templates for the reverse SSC, matching the correct template selection in
`call_ssc`.

## Steps

1. Read `kam-assemble/src/assembler.rs` in full, focusing on `call_ssc_bases`,
   `call_ssc`, and the duplex assembly call site (around line 354).

2. Add an `is_forward: bool` parameter to `call_ssc_bases`. Inside the
   function, select `template_r1` when `is_forward` and `template_r2` when
   not, mirroring the logic in `call_ssc`.

3. Update both call sites of `call_ssc_bases` to pass the correct `is_forward`
   flag.

4. Add or update a test that verifies the duplex consensus uses the correct
   template for each strand (forward SSC built from R1, reverse SSC built
   from R2).

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-001.
- The current bug means duplex consensus compares R1-vs-R1 instead of R1-vs-R2,
  reducing the error-correction benefit of duplex calling.
- Reference: `call_ssc` at lines ~429-436 shows the correct template selection.
