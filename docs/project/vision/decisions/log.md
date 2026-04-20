# Decision Log

This file is append-only. Add new entries at the bottom. Do not edit past entries.

---

## 2026-03-24: Paper framing as proof-of-concept

**Decision:** The paper presents kam as a proof-of-concept for alignment-free variant detection, not a production tool launch.

**Reasoning:** The tool is built with production quality in mind, but the paper's job is to demonstrate the approach works and is comparable. Tool launch comes later.

---

## 2026-03-24: Tumour-informed mode as primary evaluation focus

**Decision:** The paper's headline results come from tumour-informed (monitoring) mode, not discovery mode.

**Reasoning:** This is where kam achieves precision 1.0 with zero false positives. The clinical monitoring use case is the strongest story. Discovery mode results are included but secondary.

---

## 2026-03-24: config.toml as primary run interface

**Decision:** Introduce a config.toml file that configures the entire run: chemistry parameters, CLI options, thresholds, outputs. `kam run --config config.toml` becomes the primary interface.

**Reasoning:** Supports multiple chemistries without code changes. Gives users full control. Follows varforge's pattern. CLI flags still work for overrides.

---

## 2026-03-24: Expand SV types

**Decision:** Implement detection for fusions, translocations, large deletions, inversions, duplications, and tandem duplications.

**Reasoning:** The "comparable" claim requires covering the same variant types alignment-based methods handle. Current SV support is limited.

---

## 2026-03-24: Benchmark on public datasets

**Decision:** Search broadly for public datasets (not just Twist) and benchmark kam across diverse chemistries and variant types.

**Reasoning:** Demonstrates the approach generalises beyond one chemistry. Strengthens the proof-of-concept claim.

---

## 2026-03-24: Generalise beyond Twist chemistry

**Decision:** Make chemistry parameters configurable (UMI length, position, duplex/simplex, skip bases) rather than hardcoded to Twist.

**Reasoning:** The paper claims the approach works for duplex UMI sequencing generally, not just Twist. The tool should reflect that. Production use requires it.

---

## 2026-04-20: Pin bincode 1.x until v0.4.0 (RUSTSEC-2025-0141)

**Decision:** Pin `bincode = "1"` and document the advisory as allowed in `audit.toml`. Do not migrate to bincode 2.x now.

**Reasoning:** The advisory classifies bincode 1.3.3 as unmaintained, not vulnerable. All serialised data in kam is pipeline-internal — we never deserialise untrusted input, so the risk is negligible. Migrating to bincode 2.x would require an on-disk format break (or a version header and converter), which is v0.4.0 scope. Doing it now adds churn with no security benefit. Revisit at v0.4.0 planning.

---

## 2026-04-20: Accept rand 0.8 and paste 1.0.15 advisories as transitive-only (RUSTSEC-2026-0097, RUSTSEC-2024-0436)

**Decision:** Allow RUSTSEC-2026-0097 (rand 0.8 unsound) and RUSTSEC-2024-0436 (paste 1.0.15 unmaintained) in `audit.toml`. No code changes.

**Reasoning:** Both advisories are transitive — pulled in by statrs→nalgebra→simba. kam calls no rand or paste APIs directly. statrs 0.18.0 is the current latest release and still pins rand 0.8.5; there is no upgrade path available today. The unsoundness in rand 0.8 is theoretical in our context: no untrusted entropy or adversarial inputs are involved. Re-evaluate when statrs publishes a release against rand 0.9.
