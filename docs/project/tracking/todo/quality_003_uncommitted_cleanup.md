# QUALITY-003: Triage uncommitted summary TSVs and stray files

**Epic**: QUALITY
**Priority**: low
**Depends on**: none
**Status**: todo

## Goal

Decide for each uncommitted file whether it is (a) a benchmark artefact to commit, (b) a local experiment output to gitignore, or (c) junk to delete.

## Success Criteria

- [ ] `git status --short` is clean (no untracked files relevant to the project).
- [ ] `docs/benchmarking/01-snvindel/summary/` only contains TSVs referenced from a README or script.
- [ ] `.gitignore` has new entries for local-only sweep outputs if needed.
- [ ] `/update` has been run after changes.

## Steps

1. List every uncommitted path.
2. Classify each as commit / gitignore / delete.
3. Act on each classification.

## Notes

As of 2026-04-20: 18 TSVs in 01-snvindel/summary, 4 new samples_* dirs, docs/examples/explore-demo.sh, stray file `re` at project root.
