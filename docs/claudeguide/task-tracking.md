# Task Tracking

Tasks and epics live in `docs/claudetracking/`. Always check this directory
before starting any planned work.

## Folder structure

```
docs/claudetracking/
├── overallplans/   — epic-level plans (one file per major work stream)
├── todo/           — tasks ready to be picked up
├── inprogress/     — tasks currently being worked on
├── done/           — completed tasks
├── blocked/        — tasks blocked on a decision or dependency
└── open_questions.md  — design decisions needing user input
```

## Task workflow

- **Before starting**: move the task file from `todo/` to `inprogress/`.
- **After finishing**: run pre-commit checks, then move to `done/`.
- **When blocked**: move to `blocked/` and add a one-line note at the top
  stating the blocker.
- **Dependencies**: check the `Depends on` field. Do not start a task whose
  dependencies are not in `done/`.
- **Parallel work**: when tasks are independent, run them in parallel using
  separate agents. Each agent moves only its own task file.

## Task file format

```markdown
# TASK-ID: Title

**Epic**: EPIC-ID (link to overallplans file)
**Priority**: high | medium | low
**Depends on**: TASK-ID, TASK-ID (or "none")
**Status**: todo | inprogress | done | blocked

## Goal

One paragraph. What does done look like?

## Steps

1. Step one
2. Step two

## Notes

Any context, gotchas, or open decisions.
```

## Epics

Create an epic before breaking work into tasks. Never create tasks without an
epic unless the work is a one-off fix. Each epic file in `overallplans/`
lists its child tasks with IDs, file paths, and status.

## After completing work

1. Run pre-commit checks (`cargo fmt -- --check`, `cargo clippy -- -D warnings`,
   `cargo test`) and fix any failures.
2. Stage and commit with a clear message.
3. Push to the remote.

Never leave work uncommitted at the end of a session unless the user asks.

## Current epics

| Epic | File | Status |
|------|------|--------|
| BENCH-SNV | overallplans/BENCH-SNV.md | active |
| BENCH-SV | overallplans/BENCH-SV.md | active |
