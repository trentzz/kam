# FIX-001: Deduplicate parse_target_id in kam-call

**Epic**: none (standalone fix)
**Priority**: low
**Depends on**: none
**Status**: done

## Goal

Remove the duplicate `parse_target_id` function in `kam-call/src/allele.rs`
by consolidating to the `targeting.rs` version (which uses `i64` and is
needed for tolerance arithmetic). Done looks like: one `parse_target_id`
function exists in `kam-call`, exported from a shared location.

## Steps

1. Read `kam-call/src/allele.rs` and `kam-call/src/targeting.rs`.

2. Check callers of both `parse_target_id` functions. The `allele.rs` version
   returns `Option<(String, u64)>` and the `targeting.rs` version returns
   `Option<(String, i64)>`.

3. Update callers of the `allele.rs` version to use the `targeting.rs` version
   (cast `i64` to `u64` at the call site if needed, or just use `i64`).

4. Remove the `allele.rs` version.

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding NEW-006.
