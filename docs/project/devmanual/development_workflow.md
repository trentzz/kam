# Development Workflow: Claude Code Loop

## What Works Well Autonomously

**High autonomy, low risk:**
- Implementing a well-specified struct/function where interface is defined
- Writing tests for existing code
- Implementing known algorithms (quality-weighted consensus voting, Hamming distance clustering)
- Refactoring within a module without changing interfaces
- Adding structured logging and QC checks
- Writing synthetic data generator for specific chemistry once data model is defined

**Medium autonomy, needs review:**
- Implementing a new module against a defined interface
- Optimising hot paths (can find improvements, verify approach needed)
- Integrating two existing components

**Low autonomy, needs user:**
- Designing interfaces between modules
- Statistical model choices (beta-binomial vs simpler binomial)
- Any decision involving tradeoffs only user can evaluate

## Task Queue Structure

```
docs/claudeloop/
├── queue/            # Tasks to be done (picked up in filename sort order)
├── in_progress/      # Currently being worked on
├── done/             # Completed and reviewed
├── needs_review/     # Tasks Claude got stuck on
└── open_questions.md # Design questions needing user input
```

## Task File Format

Each task in `queue/` is a markdown file with precise specification:

```markdown
# Task NNN: Short Description

## Location
`kam-assemble/src/umi.rs`

## What to implement
Precise description of the function/struct/module.

## Interface (implement exactly this)
\```rust
// Exact type signatures and trait impls
\```

## Tests required
1. Test case description with expected behaviour
2. Edge case description
3. ...

## Definition of done
- `cargo test -p <crate>` passes
- `cargo clippy -p <crate> -- -D warnings` passes
- Doc comment with example
```

**The more specific the task, the more useful the output.** Vague tasks produce vague results.

## Loop Script

The loop script (`claude_loop.sh`) does:
1. Pick next task from `queue/` (sorted by filename)
2. Move to `in_progress/`
3. Run Claude Code with `--print` flag (non-interactive)
4. **Independent `cargo test` verification** (never trust self-assessment)
5. If tests pass → commit + move task to `done/`
6. If tests fail → `git stash` + move task to `needs_review/`
7. Brief pause, repeat

## Safety Boundaries

1. **Never let it touch kam-core** without an explicit task. Pre-commit hook enforces this:
   ```bash
   if git diff --cached --name-only | grep -q "kam-core/"; then
       if ! git log --format="%s" -1 | grep -q "APPROVED"; then
           echo "ERROR: Changes to kam-core require APPROVED in commit message"
           exit 1
       fi
   fi
   ```

2. **Commit after every task** — small commits for easy `git bisect`

3. **Cap tasks per session** — 5 well-specified tasks is realistic. 10+ produces hard-to-review drift.

4. **Never push automatically** — loop commits locally. Review in morning, then push.

## When Stuck Protocol

1. Design decision needed → write `// DESIGN_QUESTION: [question]` comment, move to next task
2. Tests fail after 3 attempts → write `// TODO: [what you tried]`, move on
3. Never change interface definitions in kam-core without explicit instruction
4. If unsure about scope, err on the side of smaller

## Morning Review

```bash
git log --oneline main..HEAD      # What got done
git diff main HEAD -- src/        # Quick diff review
ls docs/claudeloop/needs_review/  # Stuck tasks
cargo test                         # Full test suite confirmation
```

`needs_review/` items are morning's first tasks — understanding why it got stuck tells you if the task spec was ambiguous or if there's a genuine design question.

## Realistic Expectations

- Good night (well-specified tasks): 3-8 complete functions with tests
- Bad night (ambiguous tasks): spins on one task generating failing attempts
- 30 minutes writing precise task files > 8 hours of Claude Code with vague descriptions
