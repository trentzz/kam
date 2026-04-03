# REFACTOR-002: Deduplicate walk_paths and walk_paths_biased

**Epic**: none (standalone refactor)
**Priority**: medium
**Depends on**: none
**Status**: done

## Goal

Eliminate the code duplication between `walk_paths` and `walk_paths_biased` in
`kam-pathfind/src/walk.rs`. Both implement the same DFS with the only
difference being successor ordering. Done looks like: a shared internal
implementation parameterised by a sort function, with both public functions
calling it.

## Steps

1. Read `kam-pathfind/src/walk.rs` in full.

2. Extract a private `walk_paths_inner` function that takes a `sort_fn` closure
   (or an enum) controlling successor ordering:
   ```rust
   fn walk_paths_inner<F>(
       graph: &DeBruijnGraph,
       start: &[u8],
       end: &[u8],
       config: &WalkConfig,
       sort_successors: F,
   ) -> Vec<WalkPath>
   where
       F: Fn(&mut Vec<NodeId>),
   ```

3. Rewrite `walk_paths` and `walk_paths_biased` as thin wrappers that call
   `walk_paths_inner` with the appropriate sort closure.

4. Ensure all existing tests still pass without modification.

5. Move this task file to done and update Status.

6. Run `cargo test` and `cargo clippy -- -D warnings`. Fix all failures.

## Notes

- From codebase review finding OLD-4.1.
- Do not change the public API or behaviour — only the internal structure.
