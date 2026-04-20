# bincode 1.x Migration Strategy

**Category**: architecture
**Related epic**: DEPS

## Context

`bincode 1.3.3` is flagged UNMAINTAINED as of 2025-12-16 (RUSTSEC-2025-0141). kam uses it across:
- `kam-core/src/serialize.rs` — bincode headers + payload for all inter-stage files
- `kam-index` — HashKmerIndex persistence
- `kam-pathfind` — GraphPath + scored path artefacts
- `kam-assemble` — molecule dumps

File format compatibility is tied to the library version. `docs/project/devmanual/output_format_specs.md` documents the header format.

## Options

1. **bincode 2.x**. Breaking API change. bincode 2 is `no_std`-friendly and has an explicit `Configuration` type. All kam files written by 0.3.0 would be unreadable by a 0.4.0 that upgrades.
2. **postcard**. `no_std`, smaller output, similar API. Different format, so also breaks files.
3. **rkyv**. Zero-copy deserialisation would help kam-index loading. Large migration cost; serde-derive-style ergonomics are worse.
4. **Do nothing for now**. Advisory is "unmaintained", not "vulnerable". Pin the version and document.

## Recommendation

Option 1 (bincode 2.x) at the v0.4.0 release. Pair with a version bump in the on-disk header (`serialize.rs::Header.version`) and write a `kam convert` sub-command that reads v1 files and emits v2. Option 4 is acceptable until paper submission — the advisory is unmaintained, not unsound, and our threat model does not include untrusted serialised input (all inputs come from the same pipeline run).
