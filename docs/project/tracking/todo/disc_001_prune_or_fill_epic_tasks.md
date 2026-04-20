# DISC-2026-04-20-001: Reconcile overallplans child task lists

**Epic**: DISC-2026-04-20
**Priority**: high
**Depends on**: none
**Status**: todo

## Goal

Every epic in `overallplans/` lists child tasks with IDs and file paths. Most of those task files do not exist in `todo/`. Either create stub task files for all listed IDs or rewrite epics to list only tasks that actually exist.

## Success Criteria

- [ ] For every epic in `overallplans/`, every row in its child task table points to a file that exists in `todo/`, `inprogress/`, `done/`, or `blocked/`.
- [ ] No orphan task files (files in `todo/` not referenced by any epic).
- [ ] `/update` has been run after changes.

## Steps

1. Diff each epic's child task table against the filesystem.
2. For active epics (PAPER, PUB-BENCH, ALIGN-COMPARE, RUNTIME, CHEM-CONFIG, ML-BOOST), create stub task files with minimum success criteria.
3. For completed-but-still-listed SV sub-tasks, move to done/ or close.

## Notes

Legacy sub-tasks in `todo/` are primarily `sv_*` and `ml_*`. 21 total as of 2026-04-20.
