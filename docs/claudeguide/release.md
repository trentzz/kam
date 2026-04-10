# Release Guide

## Registry

kam publishes to **GitHub Releases** only. It does not publish to crates.io.
This is a research tool. Distribution is via compiled binary attached to the
GitHub release at https://github.com/trentzz/kam.

The crate name in Cargo.toml is `kam-bio` (to avoid collision on crates.io),
but the binary is always `kam`.

---

## Pre-release checklist

Run all checks and fix any failures before cutting the release branch.

```bash
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
```

Also verify the version number in `Cargo.toml` (`[workspace.package]`) matches
the intended release tag.

---

## Release branch

Cut from `staging` when ready to release:

```bash
git checkout -b release/vX.Y.Z staging
```

On the release branch:

1. Bump `version` in `Cargo.toml` (`[workspace.package]`).
2. Update `CHANGELOG.md` if the project maintains one.
3. Commit the version bump: `git commit -m "chore: bump version to vX.Y.Z"`.
4. Push the release branch: `git push -u origin release/vX.Y.Z`.

---

## Build the release binary

Build with the `ml` feature enabled. This is the canonical release binary and
includes the bundled ML model.

```bash
cargo build --release --features ml
```

The binary is produced at `target/release/kam`.

Strip debug symbols to reduce binary size:

```bash
strip target/release/kam
```

---

## Create the GitHub release

Use the `gh` CLI. Set `GITHUB_TOKEN` in the environment or rely on the
credentials stored by `gh auth login`.

```bash
gh release create vX.Y.Z \
  --title "kam vX.Y.Z" \
  --notes "Release notes here." \
  target/release/kam
```

If building for multiple platforms, attach each binary with a distinct name:

```bash
gh release create vX.Y.Z \
  --title "kam vX.Y.Z" \
  --notes-file RELEASE_NOTES.md \
  "target/release/kam#kam-linux-x86_64"
```

---

## Post-release steps

After the GitHub release is created:

1. Merge the release branch into `main`:

   ```bash
   git checkout main
   git merge --no-ff release/vX.Y.Z
   git push origin main
   ```

2. Tag the commit on `main`:

   ```bash
   git tag -a vX.Y.Z -m "vX.Y.Z"
   git push origin vX.Y.Z
   ```

3. Back-merge `main` into `staging` and `dev` so all branches have the version bump:

   ```bash
   git checkout staging && git merge --no-ff main && git push origin staging
   git checkout dev     && git merge --no-ff main && git push origin dev
   ```

4. Delete the release branch:

   ```bash
   git push origin --delete release/vX.Y.Z
   git branch -d release/vX.Y.Z
   ```

---

## What to exclude from source distributions

The following are excluded from source archives and should never be committed
or attached to a release:

- `bigdata/` — large generated files, benchmark outputs, model training data.
- `docs/benchmarking/` data files — scripts are included but raw data is not.
- `*.fastq.gz` — sequencing data files.
- `.env` — local environment variables.
- `target/` — build artefacts (already gitignored).

Exclusions are declared in `kam/Cargo.toml` under `[package] exclude`.

---

## Authentication

The `gh` CLI handles GitHub authentication. Before releasing, confirm you are
logged in:

```bash
gh auth status
```

If running in CI or a non-interactive environment, set `GITHUB_TOKEN` to a
personal access token with `write:packages` and `contents` scope.

---

## Version numbering

Follow semantic versioning (`MAJOR.MINOR.PATCH`).

- `PATCH` — bug fixes and minor improvements with no interface changes.
- `MINOR` — new features, new CLI flags, or new output fields. Backward compatible.
- `MAJOR` — breaking changes to the CLI, output format, or core data model.

While the project is pre-1.0, `MINOR` increments for any significant new
capability and `PATCH` for fixes.
