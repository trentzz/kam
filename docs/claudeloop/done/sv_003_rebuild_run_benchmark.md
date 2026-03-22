# Task SV-003: Rebuild and Re-Run SV Benchmark

## What to do

1. Run `cargo build --release` from the repo root.
2. Run `bash docs/benchmarking/sv/scripts/run_kam_sv.sh` to regenerate results
   for all four VAF levels (0.5%, 1%, 2%, 5%).
3. Verify Inversion calls appear in `docs/benchmarking/sv/results/kam_vaf*/variants.tsv`.

## Expected results after fix

| VAF  | DEL detected | DUP detected | INV detected |
|------|-------------|-------------|-------------|
| 5%   | ✅ PASS      | ✅ PASS      | ✅ PASS (was MNV) |
| 2%   | ✅ PASS      | ✅ PASS      | ✅ PASS or LowConfidence |
| 1%   | ✅ PASS      | ✅ PASS      | ? |
| 0.5% | ✅ PASS      | ✅ LowConfidence | ? |

## Verification

```bash
grep "Inversion" docs/benchmarking/sv/results/kam_vaf*/variants.tsv
```

Must produce at least one hit. Record the filter status and molecule count for
each VAF level.

## Definition of done

- Build succeeds.
- `grep Inversion` returns calls in at least the 5% and 2% VAF results.
- No regression: DEL and DUP still detected at all VAF levels.
