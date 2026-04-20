# DISC-2026-04-20-003: Paper scope decision

**Epic**: DISC-2026-04-20
**Priority**: high
**Depends on**: user input (see needs-review/001)
**Status**: todo

## Goal

Resolve the gap between `docs/paper/sections/05_results.tex` (which defers cross-chemistry and alignment comparison) and `docs/project/vision/paper/focus.md` (which makes both claims load-bearing).

## Success Criteria

- [ ] A dated decision entry exists in `docs/project/vision/decisions/log.md`.
- [ ] Paper sections 5.5 and 5.6 are either expanded with real results, reframed, or removed.
- [ ] Corresponding epics (PUB-BENCH, ALIGN-COMPARE) have updated priorities matching the decision.
- [ ] `/update` has been run after changes.

## Steps

1. Read `needs-review/001_defer_or_deliver_cross_chemistry.md`.
2. Discuss options with user; choose one.
3. Update paper and epic priorities accordingly.

## Notes

Blocked on user input. Recommendation in needs-review/001 is Option 3 (partial deliver): run UMI-benchmark only, parse existing thesis RaSCALL numbers into a concordance table.
