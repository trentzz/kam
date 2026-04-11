# Benchmarking Resources

This directory holds reference data files used across multiple benchmarks.

`titration.ref.v2.redux.detection.xlsx` is the alignment-based detection results for the 24-sample SNV/INDEL titration dataset (3 ng input levels: 5 ng, 15 ng, 30 ng; 8 VAF levels: 0%, 0.001%, 0.01%, 0.1%, 0.25%, 0.5%, 1%, 2%). Sheet1 contains 10,608 rows covering 375 truth variants across all non-negative-control samples. Columns are: SampleId, Chromosome, Position, Ref, Alt, Type, SampleDP, SampleAD, ifDetected(>=2supports), VAF. The `ifDetected` column flags variants where the alignment tool recorded at least 2 supporting reads (SampleAD >= 2). This file serves as the reference standard for benchmarking kam sensitivity against an alignment-based method.
