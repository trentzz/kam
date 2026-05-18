# DISC-2026-04-20-002: Make parser skip_length variable

**Epic**: DISC-2026-04-20 (also serves CHEM-CONFIG)
**Priority**: high
**Depends on**: none
**Status**: done

## Goal

`kam-assemble/src/parser.rs` currently uses `[u8; 2]` arrays for `skip_r1`/`skip_r2` and rejects any `skip_length != 2`. Generalise to any skip length (including 0).

## Success Criteria

- [x] `ParsedReadPair::skip_r1` and `skip_r2` are variable-length (Vec<u8> or similar).
- [x] `ReadStructure { umi_length: 12, skip_length: 0 }` parses without error.
- [x] A unit test for `simplex_12bp()` chemistry round-trips through the parser.
- [x] Existing Twist duplex tests still pass.
- [x] `cargo clippy --all --all-targets -- -D warnings` passes.
- [x] All tests pass.
- [x] `/update` has been run after changes.

## Steps

1. Change `skip_r1`/`skip_r2` fields in `ParsedReadPair` from `[u8; 2]` to `Vec<u8>`.
2. Update `parse_read_pair` to slice `ReadStructure::skip_length` bytes.
3. Remove the `skip_len != 2` rejection guard.
4. Add a unit test for the `simplex_12bp` preset.

## Notes

This is the remaining blocker for PUB-BENCH UMI-clustering public dataset and the CHEM-CONFIG epic.
