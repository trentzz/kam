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
