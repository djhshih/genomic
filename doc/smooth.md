# smooth port status: DNAcopy `smooth.CNA` vs C++ implementation

## Goal
Port DNAcopy's `smooth.CNA` functionality into C++, preserving the original preprocessing and smoothing semantics, and validate the result against DNAcopy on shared deterministic inputs.

## Current status
The smoothing port is implemented and tested against DNAcopy-generated expectations.

Implemented in C++ under `lib/cbs`:
- trimmed-variance threshold calculation
- DNAcopy-style finite-only preprocessing
- chromosome-local smoothing kernel
- single-profile smoothing API
- multi-sample convenience API

Public APIs:
- `cbs::smooth(...)`
- `cbs::smooth_matrix(...)`

Source files:
- `lib/cbs/smooth.hpp`
- `lib/cbs/smooth.cpp`

Reference sources used:
- `external/DNAcopy/R/DNAcopyMethods.R`
- `external/DNAcopy/R/changepoints.R`
- `external/DNAcopy/src/smoothCNA.f`

## What was ported
The port includes both layers of DNAcopy behavior.

### 1. R-side preprocessing behavior
From `smooth.CNA(...)`:
- smoothing is applied sample-by-sample
- only finite values are smoothed
- non-finite entries remain unchanged
- chromosome runs are recomputed on the finite subset only
- thresholds are derived from `trimmed.variance(...)`

### 2. Fortran smoothing kernel
From `smoothLR` in `smoothCNA.f`:
- smoothing is chromosome-local
- each point uses a radius-`k` neighborhood within its chromosome
- if any neighbor lies within `outlier.SD`, no smoothing is applied
- otherwise the local median is computed over the neighborhood
- an outlier is moved toward the neighborhood to `median ± smooth.SD`
- the replacement is not the median itself

## Current C++ API

### Single profile
```cpp
std::vector<double> cbs::smooth(const std::vector<double>& values,
                                const std::vector<int>& chrom,
                                int smooth_region = 10,
                                double outlier_sd_scale = 4.0,
                                double smooth_sd_scale = 2.0,
                                double trim = 0.025);
```

Behavior:
- `values.size() == chrom.size()` is required
- input is assumed already ordered in CNA order
- non-finite values are preserved
- chromosome labels are integer IDs

### Multi-sample convenience API
```cpp
std::vector<std::vector<double>> cbs::smooth_matrix(const std::vector<std::vector<double>>& samples,
                                                    const std::vector<int>& chrom,
                                                    int smooth_region = 10,
                                                    double outlier_sd_scale = 4.0,
                                                    double smooth_sd_scale = 2.0,
                                                    double trim = 0.025);
```

Behavior:
- each sample is smoothed independently using the same chromosome map
- this mirrors the R loop over samples in `smooth.CNA`

## Important implementation details

### Finite-only behavior
The C++ implementation mirrors DNAcopy by:
- extracting only finite values for smoothing
- rebuilding chromosome run lengths on that finite subset
- writing smoothed values back only to finite positions
- leaving `NaN`, `Inf`, and `-Inf` entries untouched

This matters because missing values change the effective local neighborhoods.

### Trimmed variance
The implementation ports DNAcopy's:
- `trimmed.variance(...)`
- `inflfact(trim)`

It uses:
- adjacent differences of the finite subset
- trimming of the sorted absolute differences
- midpoint integration for the inflation factor

### Kernel fidelity
The smoothing kernel follows the Fortran control flow closely:
- chromosome-by-chromosome processing
- neighborhood bounds clipped per chromosome
- early keep/no-smooth when any nearby point is within `outlier.SD`
- median over the full neighborhood including the focal point
- replacement to `median + smooth.SD` or `median - smooth.SD`

## Test and fixture generation
Fixture generator:
- `tests/smooth_generate.R`

Test file:
- `tests/smooth_test.cpp`

Generated fixture files:
- `tests/data/smooth_case1_input.tsv`
- `tests/data/smooth_case1_chrom.tsv`
- `tests/data/smooth_case1_expected.tsv`
- `tests/data/smooth_case2_input.tsv`
- `tests/data/smooth_case2_chrom.tsv`
- `tests/data/smooth_case2_expected.tsv`
- `tests/data/smooth_case3_input.tsv`
- `tests/data/smooth_case3_chrom.tsv`
- `tests/data/smooth_case3_expected.tsv`

These fixtures are generated from actual DNAcopy `smooth.CNA(...)` calls.

## What the tests cover
`tests/smooth_test.cpp` validates full output equality against DNAcopy-generated expectations.

Coverage includes:
- single-profile smoothing
- missing/non-finite handling
- chromosome-boundary behavior
- multi-sample smoothing

Comparisons are made on the complete output vectors, not spot checks.

## Build and test integration
Integrated into the project build via:
- `CMakeLists.txt`
- `tests/CMakeLists.txt`

The test generator is run from CMake via:
- `Rscript tests/smooth_generate.R`

The regression test target is:
- `smooth_test`

## Current conclusion
The C++ smoothing port should now be described as matching DNAcopy `smooth.CNA` for the regression cases currently covered by the test suite.

What is true now:
- preprocessing semantics are matched
- Fortran smoothing semantics are matched
- DNAcopy-generated fixtures are used for regression
- non-finite and chromosome-boundary behavior are covered
- single-profile and multi-sample behavior are covered
- full `ctest` passes

## Remaining work
No active blocker remains for the currently covered smoothing scope.

Possible future work:
- broaden the fixture matrix with more edge cases
- add stress tests for short chromosomes and tiny finite subsets
- add higher-level integration if a native C++ CNA object abstraction is introduced

## Bottom line
The `smooth.CNA` port is implemented, validated against DNAcopy, and integrated into the library and test suite.
