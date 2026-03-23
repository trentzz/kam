# FIX-003: Fix infer_k_from_entries stub and add --kmer-size to pathfind

**Epic**: none (standalone fix)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Fix two related issues in `kam/src/commands/pathfind.rs`:
1. `infer_k_from_entries` always returns `Some(31)` regardless of the index
   content (stub implementation).
2. The `pathfind` subcommand has no `--kmer-size` flag, so there is no way
   to override the inferred k.

Done looks like: `infer_k_from_entries` reads the k-mer size from the index
entries (or returns `None` if the index is empty), and `PathfindArgs` has a
`--kmer-size` flag that overrides inference.

## Steps

1. Read `kam/src/commands/pathfind.rs` and `kam/src/cli.rs` (PathfindArgs).

2. Read the index file format (check how k-mer size is stored in the bincode
   index — look at `kam-index/src/` and the FileType::KmerIndex entries).

3. Fix `infer_k_from_entries`: iterate entries to find the actual kmer length
   from the stored keys, or read it from the file header if available. If the
   index is empty, return `None`.

4. Add `--kmer-size: Option<usize>` to `PathfindArgs` in `kam/src/cli.rs`.
   When provided, use it directly; otherwise fall back to `infer_k_from_entries`.

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review findings OLD-1.4 and OLD-8.2.
- The `run` subcommand passes k explicitly and is not affected.
- Only the standalone `kam pathfind` subcommand is broken.
