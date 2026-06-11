# CBS port status: DNAcopy vs C++ implementation

## Goal
Port DNAcopy's CBS implementation into C++, including weighted CBS, preserve the original algorithmic behavior, and validate the result against DNAcopy on shared deterministic inputs.

## Current status
The CBS port is now implemented and tested end-to-end against DNAcopy for the regression cases covered in `tests/cbs_test.cpp` with respect to the root of the project.

Implemented in C++ under `lib/cbs`:
- low-level unweighted CBS routines
- low-level weighted CBS routines
- top-level segmentation drivers
- pruning logic needed to match DNAcopy final segmentation behavior

Public APIs now include:
- `cbs::SegmentationResult`
- `cbs::segment(...)`
- `cbs::segment_weighted(...)`

## Shared input and expected-output generation
`tests/cbs_compare.R` is the single CBS fixture generator.

It writes:
- raw input profiles
- optional weights
- expected DNAcopy segmentations

Generated files include:
- `tests/data/cbs_case1_input.tsv`
- `tests/data/cbs_case1_expected.tsv`
- `tests/data/cbs_case1_perm_alt_expected.tsv`
- `tests/data/cbs_case1_hybrid_expected.tsv`
- `tests/data/cbs_case1_hybrid_alt_expected.tsv`
- `tests/data/cbs_case2_weighted_input.tsv`
- `tests/data/cbs_case2_weighted_weights.tsv`
- `tests/data/cbs_case2_weighted_expected.tsv`
- `tests/data/cbs_case2_weighted_perm_alt_expected.tsv`
- `tests/data/cbs_case2_weighted_hybrid_expected.tsv`
- `tests/data/cbs_case2_weighted_hybrid_alt_expected.tsv`
- `tests/data/cbs_case3_noisy_input.tsv`
- `tests/data/cbs_case3_noisy_expected.tsv`
- `tests/data/cbs_case3_noisy_perm_alt_expected.tsv`
- `tests/data/cbs_case3_noisy_hybrid_expected.tsv`
- `tests/data/cbs_case3_noisy_hybrid_alt_expected.tsv`
- `tests/data/cbs_case4_noisy_input.tsv`
- `tests/data/cbs_case4_noisy_expected.tsv`
- `tests/data/cbs_case4_noisy_perm_alt_expected.tsv`
- `tests/data/cbs_case4_noisy_hybrid_expected.tsv`
- `tests/data/cbs_case4_noisy_hybrid_alt_expected.tsv`

The noisy examples are reproducible via fixed seeds in `tests/cbs_compare.R`.

The script supports:
- fixture generation: `Rscript tests/cbs_compare.R`
- hybrid inspection output: `Rscript tests/cbs_compare.R print-hybrid`

## What the C++ tests now check
`tests/cbs_test.cpp` reads the generated TSV files directly, so DNAcopy expectations and the C++ implementation are compared on exactly the same inputs.

The tests intentionally separate two semantic layers.

### 1. Raw-kernel behavior
These tests validate literal low-level behavior against the original Fortran-style primitives:
- `tmaxo(...)`
- `wtmaxo(...)`

This is important because the maximizing interval returned by the raw statistic kernel is not the same thing as a final segmented profile.

### 2. Final segmentation-driver behavior
These tests validate:
- `cbs::segment(...)`
- `cbs::segment_weighted(...)`

against DNAcopy expected outputs, comparing:
- full segment-length vectors
- full segment-mean vectors

## Parameter coverage currently exercised
The regression suite covers multiple DNAcopy-derived parameterizations.

### Simple unweighted case
Checked for:
- permutation mode, default settings
- permutation mode, alternate settings
- hybrid mode, default settings
- hybrid mode, alternate settings

### Simple weighted case
Checked for:
- permutation mode, default settings
- permutation mode, alternate settings
- hybrid mode, default settings
- hybrid mode, alternate settings

### Noisy unweighted cases
Checked for both noisy profiles:
- permutation mode, default settings
- permutation mode, alternate settings
- hybrid mode, default settings
- hybrid mode, alternate settings

## Important implementation findings and fixes
Several correctness issues were identified and fixed during the port.

### Critical permutation fast-accept bug
In the C++ port, `tpermp(...)` and `wtpermp(...)` initially returned `1.0` instead of `0.0` in the strong-signal fast-accept branch.

That was wrong relative to DNAcopy / Fortran semantics and caused real changepoints to be rejected.

This is now fixed.

### Weighted cumulative-weight construction
The weighted driver now uses the DNAcopy-compatible construction:
- `cwts = cumsum(weights) / sqrt(sum(weights))`

### Permutation buffer reset
Permutation code now reinitializes the permuted working buffer on each iteration, matching the intended behavior.

### Final changepoint path alignment
`fndcpt(...)` / `wfindcpt(...)` and the top-level segmentation path were corrected so the final segmentation behavior matches DNAcopy rather than only matching raw maximizing-interval kernels.

## Current conclusion
The port should now be described as matching DNAcopy end-to-end for the regression cases currently covered by the test suite.

What is true now:
- shared DNAcopy/C++ inputs are used
- expected outputs come directly from DNAcopy
- low-level kernel behavior is tested separately from final segmentation behavior
- unweighted and weighted drivers are implemented
- noisy and simple profiles are covered
- hybrid and permutation modes are covered
- full `ctest` passes

## Remaining work
No active blocker remains for CBS in the currently covered scope.

Possible future work:
- broaden the regression matrix further
- add more weighted noisy fixtures
- add more parameter combinations if desired
- keep fixture provenance synchronized with `tests/cbs_compare.R`

## Bottom line
This project is past the earlier gap-analysis stage.

The original architectural gap has been closed:
- the C++ port now includes a top-level CBS segmentation driver
- the regression tests validate full DNAcopy-equivalent behavior on the current fixture set
- `tests/cbs_compare.R` is the single source of fixture generation and hybrid inspection
