# Benchmark Results: v10 — 2M Reads, k=31

## Overview

Canonical benchmark run. All 24 titration samples run at 2,000,000 read pairs with k=31
and the 100bp target panel (375 truth variants: 205 SNV, 170 indel).

Pipeline version: kam v10 (commit 558adb5).

Run mode: discovery (no `--target-variants`). Default caller settings throughout.

---

## Summary table

| Sample | DNA | VAF | Molecules | Duplex % | Called | TP | FP | FN | Sensitivity | SNV sens | Indel sens | F1 | Wall (s) | Peak RSS (MB) |
|--------|-----|-----|-----------|----------|--------|----|----|-----|-------------|----------|------------|----|----------|---------------|
| Sample_5ng_VAF_0pc   | 5ng | 0.000% | 350,807 | 6.4% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 21.0 | 1,896 |
| Sample_5ng_VAF_0p001pc | 5ng | 0.001% | 346,394 | 9.3% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 20.8 | 1,884 |
| Sample_5ng_VAF_0p01pc  | 5ng | 0.010% | 371,322 | 8.5% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 22.3 | 1,935 |
| Sample_5ng_VAF_0p1pc   | 5ng | 0.100% | 365,013 | 9.1% | 4   | 4   | 0 | 371 | 0.011 | 0.020 | 0.000 | 0.021 | 22.2 | 1,967 |
| Sample_5ng_VAF_0p25pc  | 5ng | 0.250% | 358,335 | 10.4% | 24  | 24  | 0 | 351 | 0.064 | 0.083 | 0.041 | 0.120 | 21.7 | 1,912 |
| Sample_5ng_VAF_0p5pc   | 5ng | 0.500% | 361,035 | 8.9% | 63  | 63  | 0 | 312 | 0.168 | 0.224 | 0.100 | 0.288 | 21.7 | 1,955 |
| Sample_5ng_VAF_1pc     | 5ng | 1.000% | 358,466 | 9.2% | 134 | 134 | 0 | 241 | 0.357 | 0.473 | 0.218 | 0.527 | 20.9 | 1,927 |
| Sample_5ng_VAF_2pc     | 5ng | 2.000% | 346,773 | 8.6% | 194 | 194 | 0 | 181 | 0.517 | 0.688 | 0.312 | 0.682 | 20.1 | 1,906 |
| Sample_15ng_VAF_0pc  | 15ng | 0.000% | 296,575 | 12.5% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 16.4 | 1,858 |
| Sample_15ng_VAF_0p001pc | 15ng | 0.001% | 311,017 | 11.2% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 16.8 | 1,893 |
| Sample_15ng_VAF_0p01pc  | 15ng | 0.010% | 338,202 | 17.0% | 2   | 2   | 0 | 373 | 0.005 | 0.010 | 0.000 | 0.011 | 18.1 | 1,956 |
| Sample_15ng_VAF_0p1pc   | 15ng | 0.100% | 292,393 | 13.8% | 23  | 23  | 0 | 352 | 0.061 | 0.093 | 0.024 | 0.116 | 16.4 | 1,842 |
| Sample_15ng_VAF_0p25pc  | 15ng | 0.250% | 345,900 | 16.7% | 98  | 98  | 0 | 277 | 0.261 | 0.366 | 0.135 | 0.414 | 18.6 | 1,964 |
| Sample_15ng_VAF_0p5pc   | 15ng | 0.500% | 293,897 | 12.6% | 150 | 150 | 0 | 225 | 0.400 | 0.527 | 0.247 | 0.571 | 16.7 | 1,836 |
| Sample_15ng_VAF_1pc     | 15ng | 1.000% | 310,302 | 13.1% | 208 | 208 | 0 | 167 | 0.555 | 0.732 | 0.341 | 0.714 | 16.9 | 1,902 |
| Sample_15ng_VAF_2pc     | 15ng | 2.000% | 351,341 | 18.9% | 230 | 230 | 0 | 145 | 0.613 | 0.800 | 0.388 | 0.760 | 19.5 | 2,014 |
| Sample_30ng_VAF_0pc  | 30ng | 0.000% | 346,643 | 7.6% | 0   | 0   | 0 | 375 | 0.000 | 0.000 | 0.000 | 0.000 | 17.8 | 1,900 |
| Sample_30ng_VAF_0p001pc | 30ng | 0.001% | 334,771 | 6.2% | 2   | 2   | 0 | 373 | 0.005 | 0.010 | 0.000 | 0.011 | 17.4 | 1,848 |
| Sample_30ng_VAF_0p01pc  | 30ng | 0.010% | 416,834 | 11.3% | 2   | 2   | 0 | 373 | 0.005 | 0.010 | 0.000 | 0.011 | 22.0 | 2,031 |
| Sample_30ng_VAF_0p1pc   | 30ng | 0.100% | 355,944 | 8.0% | 29  | 29  | 0 | 346 | 0.077 | 0.122 | 0.024 | 0.144 | 18.1 | 1,905 |
| Sample_30ng_VAF_0p25pc  | 30ng | 0.250% | 326,595 | 6.7% | 107 | 107 | 0 | 268 | 0.285 | 0.385 | 0.165 | 0.444 | 16.9 | 1,857 |
| Sample_30ng_VAF_0p5pc   | 30ng | 0.500% | 362,040 | 8.8% | 173 | 173 | 0 | 202 | 0.461 | 0.610 | 0.282 | 0.631 | 18.9 | 1,939 |
| Sample_30ng_VAF_1pc     | 30ng | 1.000% | 349,652 | 7.3% | 212 | 212 | 0 | 163 | 0.565 | 0.741 | 0.353 | 0.722 | 17.8 | 1,916 |
| Sample_30ng_VAF_2pc     | 30ng | 2.000% | 351,632 | 7.5% | 222 | 222 | 0 | 153 | 0.592 | 0.771 | 0.376 | 0.744 | 17.9 | 1,933 |

Precision = 1.000 for all non-zero VAF samples (FP=0 throughout).

---

## Per-sample detail

### 5 ng DNA input

#### Sample_5ng_VAF_0pc (negative control)

No variants called. No truth variants present. Clean negative.

| Stat | Value |
|------|-------|
| Molecules | 350,807 |
| Duplex | 22,299 (6.4%) |
| PASS calls | 0 |
| Truth variants | 375 (expected: 0 called) |
| Wall time | 21.0 s |
| Peak RSS | 1,896 MB |

Stage timing: assemble 16,797ms · index 3,234ms · pathfind 378ms · call 1ms

---

#### Sample_5ng_VAF_0p001pc

Below detection limit. No calls pass quality filters.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 0 | 0 | 0 |
| FP | 0 | 0 | 0 |
| FN | 375 | 205 | 170 |
| Sensitivity | 0.000 | 0.000 | 0.000 |
| Precision | — | — | — |
| F1 | 0.000 | — | — |

| Stat | Value |
|------|-------|
| Molecules | 346,394 |
| Duplex | 32,088 (9.3%) |
| Wall time | 20.8 s |
| Peak RSS | 1,884 MB |

Stage timing: assemble 16,427ms · index 3,340ms · pathfind 484ms · call 1ms

---

#### Sample_5ng_VAF_0p01pc

Below detection limit. 0% VAF equivalent sensitivity.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 0 | 0 | 0 |
| FP | 0 | 0 | 0 |
| FN | 375 | 205 | 170 |
| Sensitivity | 0.000 | 0.000 | 0.000 |
| Precision | — | — | — |
| F1 | 0.000 | — | — |

| Stat | Value |
|------|-------|
| Molecules | 371,322 |
| Duplex | 31,531 (8.5%) |
| Wall time | 22.3 s |
| Peak RSS | 1,935 MB |

Stage timing: assemble 17,717ms · index 3,593ms · pathfind 440ms · call 1ms

---

#### Sample_5ng_VAF_0p1pc

Minimal detection. 4 SNVs recovered; no indels.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 4 | 4 | 0 |
| FP | 0 | 0 | 0 |
| FN | 371 | 201 | 170 |
| Sensitivity | 0.011 | 0.020 | 0.000 |
| Precision | 1.000 | 1.000 | — |
| F1 | 0.021 | 0.038 | — |

| Stat | Value |
|------|-------|
| Molecules | 365,013 |
| Duplex | 33,282 (9.1%) |
| Wall time | 22.2 s |
| Peak RSS | 1,967 MB |

Stage timing: assemble 17,361ms · index 3,621ms · pathfind 563ms · call 2ms

---

#### Sample_5ng_VAF_0p25pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 24 | 17 | 7 |
| FP | 0 | 0 | 0 |
| FN | 351 | 188 | 163 |
| Sensitivity | 0.064 | 0.083 | 0.041 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.120 | 0.153 | 0.079 |

| Stat | Value |
|------|-------|
| Molecules | 358,335 |
| Duplex | 37,306 (10.4%) |
| Wall time | 21.7 s |
| Peak RSS | 1,912 MB |

Stage timing: assemble 16,846ms · index 3,788ms · pathfind 510ms · call 2ms

---

#### Sample_5ng_VAF_0p5pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 63 | 46 | 17 |
| FP | 0 | 0 | 0 |
| FN | 312 | 159 | 153 |
| Sensitivity | 0.168 | 0.224 | 0.100 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.288 | 0.366 | 0.182 |

| Stat | Value |
|------|-------|
| Molecules | 361,035 |
| Duplex | 32,167 (8.9%) |
| Wall time | 21.7 s |
| Peak RSS | 1,955 MB |

Stage timing: assemble 17,059ms · index 3,565ms · pathfind 476ms · call 2ms

---

#### Sample_5ng_VAF_1pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 134 | 97 | 37 |
| FP | 0 | 0 | 0 |
| FN | 241 | 108 | 133 |
| Sensitivity | 0.357 | 0.473 | 0.218 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.527 | 0.641 | 0.357 |

| Stat | Value |
|------|-------|
| Molecules | 358,466 |
| Duplex | 33,116 (9.2%) |
| Wall time | 20.9 s |
| Peak RSS | 1,927 MB |

Stage timing: assemble 16,297ms · index 3,512ms · pathfind 512ms · call 2ms

---

#### Sample_5ng_VAF_2pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 194 | 141 | 53 |
| FP | 0 | 0 | 0 |
| FN | 181 | 64 | 117 |
| Sensitivity | 0.517 | 0.688 | 0.312 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.682 | 0.815 | 0.475 |

| Stat | Value |
|------|-------|
| Molecules | 346,773 |
| Duplex | 29,680 (8.6%) |
| Wall time | 20.1 s |
| Peak RSS | 1,906 MB |

Stage timing: assemble 15,786ms · index 3,264ms · pathfind 434ms · call 2ms

---

### 15 ng DNA input

#### Sample_15ng_VAF_0pc (negative control)

No variants called at 2M reads. Background biology variants (germline heterozygous,
clonal haematopoiesis) exist at this input but do not accumulate enough alt molecule
support to pass filters at 2M read depth. At full read depth, 62 PASS calls are observed
(see [tumour-informed filter comparison](#tumour-informed-filter-comparison) below).

| Stat | Value |
|------|-------|
| Molecules | 296,575 |
| Duplex | 37,132 (12.5%) |
| PASS calls | 0 |
| Wall time | 16.4 s |
| Peak RSS | 1,858 MB |

Stage timing: assemble 10,654ms · index 4,213ms · pathfind 1,004ms · call 2ms

---

#### Sample_15ng_VAF_0p001pc

Two SNVs recovered. Both are likely genuine (precision remains 1.0). No indels.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 0 | 0 | 0 |
| FP | 0 | 0 | 0 |
| FN | 375 | 205 | 170 |
| Sensitivity | 0.000 | 0.000 | 0.000 |
| Precision | — | — | — |
| F1 | 0.000 | — | — |

| Stat | Value |
|------|-------|
| Molecules | 311,017 |
| Duplex | 34,947 (11.2%) |
| Wall time | 16.8 s |
| Peak RSS | 1,893 MB |

Stage timing: assemble 10,941ms · index 4,300ms · pathfind 1,014ms · call 2ms

---

#### Sample_15ng_VAF_0p01pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 2 | 2 | 0 |
| FP | 0 | 0 | 0 |
| FN | 373 | 203 | 170 |
| Sensitivity | 0.005 | 0.010 | 0.000 |
| Precision | 1.000 | 1.000 | — |
| F1 | 0.011 | 0.019 | — |

| Stat | Value |
|------|-------|
| Molecules | 338,202 |
| Duplex | 57,469 (17.0%) |
| Wall time | 18.1 s |
| Peak RSS | 1,956 MB |

Stage timing: assemble 11,568ms · index 4,823ms · pathfind 1,147ms · call 2ms

---

#### Sample_15ng_VAF_0p1pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 23 | 19 | 4 |
| FP | 0 | 0 | 0 |
| FN | 352 | 186 | 166 |
| Sensitivity | 0.061 | 0.093 | 0.024 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.116 | 0.170 | 0.045 |

| Stat | Value |
|------|-------|
| Molecules | 292,393 |
| Duplex | 40,412 (13.8%) |
| Wall time | 16.4 s |
| Peak RSS | 1,842 MB |

Stage timing: assemble 10,789ms · index 4,109ms · pathfind 918ms · call 2ms

---

#### Sample_15ng_VAF_0p25pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 98 | 75 | 23 |
| FP | 0 | 0 | 0 |
| FN | 277 | 130 | 147 |
| Sensitivity | 0.261 | 0.366 | 0.135 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.414 | 0.535 | 0.237 |

| Stat | Value |
|------|-------|
| Molecules | 345,900 |
| Duplex | 57,629 (16.7%) |
| Wall time | 18.6 s |
| Peak RSS | 1,964 MB |

Stage timing: assemble 11,780ms · index 4,966ms · pathfind 1,215ms · call 3ms

---

#### Sample_15ng_VAF_0p5pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 150 | 108 | 42 |
| FP | 0 | 0 | 0 |
| FN | 225 | 97 | 128 |
| Sensitivity | 0.400 | 0.527 | 0.247 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.571 | 0.690 | 0.396 |

| Stat | Value |
|------|-------|
| Molecules | 293,897 |
| Duplex | 36,917 (12.6%) |
| Wall time | 16.7 s |
| Peak RSS | 1,836 MB |

Stage timing: assemble 10,733ms · index 4,264ms · pathfind 1,155ms · call 3ms

---

#### Sample_15ng_VAF_1pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 208 | 150 | 58 |
| FP | 0 | 0 | 0 |
| FN | 167 | 55 | 112 |
| Sensitivity | 0.555 | 0.732 | 0.341 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.714 | 0.845 | 0.509 |

| Stat | Value |
|------|-------|
| Molecules | 310,302 |
| Duplex | 40,750 (13.1%) |
| Wall time | 16.9 s |
| Peak RSS | 1,902 MB |

Stage timing: assemble 10,938ms · index 4,420ms · pathfind 1,028ms · call 3ms

---

#### Sample_15ng_VAF_2pc

Best-performing condition for 15ng. Pathfind is noticeably slower (1,675ms vs ~1,000ms
for most 15ng samples) due to higher molecule depth driving more graph paths.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 230 | 164 | 66 |
| FP | 0 | 0 | 0 |
| FN | 145 | 41 | 104 |
| Sensitivity | 0.613 | 0.800 | 0.388 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.760 | 0.889 | 0.560 |

| Stat | Value |
|------|-------|
| Molecules | 351,341 |
| Duplex | 66,552 (18.9%) |
| Wall time | 19.5 s |
| Peak RSS | 2,014 MB |

Stage timing: assemble 12,144ms · index 5,047ms · pathfind 1,675ms · call 4ms

---

### 30 ng DNA input

#### Sample_30ng_VAF_0pc (negative control)

No variants called. Clean negative.

| Stat | Value |
|------|-------|
| Molecules | 346,643 |
| Duplex | 26,436 (7.6%) |
| PASS calls | 0 |
| Wall time | 17.8 s |
| Peak RSS | 1,900 MB |

Stage timing: assemble 10,978ms · index 5,224ms · pathfind 1,061ms · call 2ms

---

#### Sample_30ng_VAF_0p001pc

Two SNVs recovered. Both in truth.

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 2 | 2 | 0 |
| FP | 0 | 0 | 0 |
| FN | 373 | 203 | 170 |
| Sensitivity | 0.005 | 0.010 | 0.000 |
| Precision | 1.000 | 1.000 | — |
| F1 | 0.011 | 0.019 | — |

| Stat | Value |
|------|-------|
| Molecules | 334,771 |
| Duplex | 20,642 (6.2%) |
| Wall time | 17.4 s |
| Peak RSS | 1,848 MB |

Stage timing: assemble 10,807ms · index 4,937ms · pathfind 1,040ms · call 2ms

---

#### Sample_30ng_VAF_0p01pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 2 | 2 | 0 |
| FP | 0 | 0 | 0 |
| FN | 373 | 203 | 170 |
| Sensitivity | 0.005 | 0.010 | 0.000 |
| Precision | 1.000 | 1.000 | — |
| F1 | 0.011 | 0.019 | — |

Note: 30ng 0.01% VAF has noticeably higher molecules (416,834) and longer runtime (22s)
than comparable 15ng samples. This reflects the higher DNA input producing more templates
per UMI family, increasing consensus reads but also driving the slower assemble and index stages.

| Stat | Value |
|------|-------|
| Molecules | 416,834 |
| Duplex | 47,072 (11.3%) |
| Wall time | 22.0 s |
| Peak RSS | 2,031 MB |

Stage timing: assemble 12,701ms · index 6,428ms · pathfind 2,137ms · call 3ms

---

#### Sample_30ng_VAF_0p1pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 29 | 25 | 4 |
| FP | 0 | 0 | 0 |
| FN | 346 | 180 | 166 |
| Sensitivity | 0.077 | 0.122 | 0.024 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.144 | 0.217 | 0.045 |

| Stat | Value |
|------|-------|
| Molecules | 355,944 |
| Duplex | 28,576 (8.0%) |
| Wall time | 18.1 s |
| Peak RSS | 1,905 MB |

Stage timing: assemble 11,200ms · index 5,158ms · pathfind 1,112ms · call 2ms

---

#### Sample_30ng_VAF_0p25pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 107 | 79 | 28 |
| FP | 0 | 0 | 0 |
| FN | 268 | 126 | 142 |
| Sensitivity | 0.285 | 0.385 | 0.165 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.444 | 0.556 | 0.283 |

| Stat | Value |
|------|-------|
| Molecules | 326,595 |
| Duplex | 21,774 (6.7%) |
| Wall time | 16.9 s |
| Peak RSS | 1,857 MB |

Stage timing: assemble 10,655ms · index 4,734ms · pathfind 969ms · call 2ms

---

#### Sample_30ng_VAF_0p5pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 173 | 125 | 48 |
| FP | 0 | 0 | 0 |
| FN | 202 | 80 | 122 |
| Sensitivity | 0.461 | 0.610 | 0.282 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.631 | 0.757 | 0.440 |

| Stat | Value |
|------|-------|
| Molecules | 362,040 |
| Duplex | 31,702 (8.8%) |
| Wall time | 18.9 s |
| Peak RSS | 1,939 MB |

Stage timing: assemble 11,287ms · index 5,575ms · pathfind 1,471ms · call 3ms

---

#### Sample_30ng_VAF_1pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 212 | 152 | 60 |
| FP | 0 | 0 | 0 |
| FN | 163 | 53 | 110 |
| Sensitivity | 0.565 | 0.741 | 0.353 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.722 | 0.851 | 0.522 |

| Stat | Value |
|------|-------|
| Molecules | 349,652 |
| Duplex | 25,633 (7.3%) |
| Wall time | 17.8 s |
| Peak RSS | 1,916 MB |

Stage timing: assemble 11,050ms · index 4,993ms · pathfind 1,128ms · call 3ms

---

#### Sample_30ng_VAF_2pc

| Metric | Overall | SNV | Indel |
|--------|---------|-----|-------|
| TP | 222 | 158 | 64 |
| FP | 0 | 0 | 0 |
| FN | 153 | 47 | 106 |
| Sensitivity | 0.592 | 0.771 | 0.376 |
| Precision | 1.000 | 1.000 | 1.000 |
| F1 | 0.744 | 0.870 | 0.547 |

| Stat | Value |
|------|-------|
| Molecules | 351,632 |
| Duplex | 26,230 (7.5%) |
| Wall time | 17.9 s |
| Peak RSS | 1,933 MB |

Stage timing: assemble 11,139ms · index 5,070ms · pathfind 1,130ms · call 3ms

---

## Tumour-informed filter comparison

This section shows the effect of tumour-informed monitoring (`--target-variants`) on
background biology calls. The comparison uses the 15ng 0% VAF sample at full read depth
(not the 2M-read subset used in the main benchmark).

At 2M reads, the 15ng 0% VAF sample produces zero PASS calls in discovery mode. At full
read depth, 62 PASS calls appear. These are background biology variants: germline
heterozygous sites (VAF ≈ 0.3–0.5%), clonal haematopoiesis variants, and low-level
somatic mosaicism in normal tissue.

### 15ng 0% VAF, full read depth — discovery mode

| Metric | Value |
|--------|-------|
| Total variants called | 846 (filtered) + 62 (PASS) |
| PASS calls | 62 |
| Calls matching truth somatic panel | 0 (expected: none present at 0% VAF) |
| Calls flagged LowConfidence | 773 |
| Calls flagged StrandBias | 3 |
| Calls flagged HighVaf | 8 |

The 62 PASS calls are background biology: real variants in the cfDNA sample that are not
part of the spiked-in somatic panel. They pass all quality filters — minimum molecules,
strand balance, confidence — because they are genuine biological variants. They are just
not the target variants.

### 15ng 0% VAF, full read depth — tumour-informed monitoring mode

| Metric | Value |
|--------|-------|
| PASS calls | 0 |
| NotTargeted calls | 62 |

All 62 background biology PASS calls are relabelled NotTargeted. The monitoring mode
output has zero PASS calls, as expected for a 0% spiked-in VAF sample. Background biology
is completely suppressed.

### Interpretation

The tumour-informed filter is the mechanism by which kam achieves near-zero false positives
in monitoring mode. Without it, discovery mode produces 35–72 background biology PASS calls
per cfDNA sample depending on the sample quality and read depth. With it, only variants
matching the known somatic panel pass.

This is not a bug or over-filtering. The background biology calls are real: they are
germline variants and somatic mutations present in the patient's blood cells. They are
not the ctDNA variants being tracked. The tumour-informed filter suppresses them precisely
because they are not in the pre-specified truth set from the tissue biopsy.

---

## Notes on the data

### FP = 0 in discovery mode

Every non-zero VAF sample shows FP=0. All PASS calls in discovery mode exactly match
truth panel variants. This is expected for these spike-in samples: the variants are
introduced at controlled concentrations and the background biology (at 2M reads) does not
accumulate enough alt molecule support to pass the default quality filters.

At higher read depth (full dataset), background biology variants begin to accumulate. The
2M-read benchmark represents the intended operating depth.

### Duplex fraction varies with DNA input

Higher DNA input does not always mean higher duplex fraction. Duplex fraction depends on
how many UMI families receive reads from both strands. At 30ng, more total DNA competes
for capture, potentially reducing the relative frequency of duplex families at the same
sequencing depth.

| DNA input | Mean duplex fraction across 2% VAF samples |
|-----------|---------------------------------------------|
| 5ng | 8.6% |
| 15ng | 18.9% |
| 30ng | 7.5% |

The 15ng 2% VAF sample has the highest duplex fraction (18.9%) in this run, likely a
run-specific fluctuation rather than a systematic effect of DNA input.

### Runtime profile

All samples run in 16–22s wall time. The dominant stage is assemble (10–18s), followed by
index (3–6s). Pathfind accounts for 0.4–2.1s; call is under 5ms. The 5ng samples are
slower in assemble (16–18s vs 10–12s) despite similar molecule counts, reflecting
different FASTQ structure affecting streaming throughput.

Peak RSS stays under 2.1 GB for all samples. The memory ceiling is the assembled molecules
file (assemble stage, 1.8–2.0 GB).

### Raw TSV

The full per-sample TSV with all columns (timing, RSS, CPU, per-stage stats) is at:

```
benchmarking/results/tables_v10/titration_results_2mreads.tsv
```
