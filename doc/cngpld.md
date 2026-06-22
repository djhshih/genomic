# cngpld summarize_cn port

This project includes a C++ port of the `summarize_cn` logic from:

- `external/cngpld/R/seg.R`

## Scope

Implemented:
- `summarize_cn_at_position`
- `summarize_cn`

Not implemented here:
- full GRanges support
- `compare_segs`
- `count_cn`
- `count_segs`
- centering helpers
- chromosome-arm helpers
- GPLDIFF fitting pipeline

## Files

Library:
- `lib/cngpld/summarize.hpp`
- `lib/cngpld/summarize.cpp`

Tests:
- `tests/cngpld_generate.R`
- `tests/cngpld_test.cpp`

## Data model

The port reuses existing repository classes rather than introducing a duplicate segment record type.

Input is provided as:
- `cna::SegmentedSampleSet<rvalue>`
- selected by sample name and chromosome index

Minimal read-only accessor support added for this:
- `const SegmentedSample* cna::SegmentedSampleSet<V>::sample(const std::string&) const`
- const chromosome accessors on `cna::Sample<Chromosome>`

No `getSamples()` accessor was added.
Traversal continues to use the repository’s existing iterator style.

## API

```cpp
namespace cngpld {

struct CNSummaryPoint {
    position pos;
    double value;
};

using CNSummary = std::vector<CNSummaryPoint>;

double summarize_cn_at_position(
    const cna::SegmentedSampleSet<rvalue>& seg,
    const std::string& sample,
    size_t chrom_index,
    position pos,
    int direction,
    double cutoff);

CNSummary summarize_cn(
    const cna::SegmentedSampleSet<rvalue>& seg,
    const std::string& sample,
    size_t chrom_index,
    int direction,
    double cutoff,
    const std::vector<position>* positions = NULL);

}
```

## Semantics

The implementation follows `external/cngpld/R/seg.R`.

### Default positions
If `positions == NULL`, the queried positions are:
- all segment starts
- all segment ends
- sorted and uniqued

Equivalent R logic:

```r
sort(unique(c(start(gr), end(gr))))
```

### Overlap rule
A segment overlaps a queried position if:

```text
segment.start <= pos && pos <= segment.end
```

This is inclusive at both ends.

### Direction and thresholding
For each overlapping segment, compute:

```text
affected = direction * seg.value
```

Where:
- `direction = 1` for amplification summary
- `direction = -1` for deletion summary

A segment contributes to the numerator only if:

```text
affected > cutoff
```

Thresholding is strict.
Equality does not pass.

### Returned value at a position
If at least one overlapping segment passes threshold:

```text
sum(exp(affected for passing overlaps)) / number_of_all_overlaps
```

Otherwise:

```text
0
```

This denominator rule is important:
- denominator = all overlapping segments
- not only altered overlaps

That matches the original R logic exactly.

## Validation

Regression fixtures are generated in R by `tests/cngpld_generate.R`.

Because the original R code path depends on Bioconductor GRanges APIs that are not available in this environment, the fixture generator uses a direct R reimplementation of the exact `summarize_cn` / `summarize_cn_at_position` logic from `seg.R`.

Covered test cases:
- amplification path
- deletion path
- denominator semantics
- explicit queried positions
- empty-overlap positions returning `0`

## Notes

- chromosome selection is explicit via `chrom_index`
- sample selection is explicit via sample name
- invalid `direction` throws
- unknown sample throws
- malformed segments with `start > end` throw
