# TODO

## Conversion features

- [ ] Implement segmented-to-raw conversion in `lib/RawSampleSet.hpp`
  - `RawSampleSet<V>::RawSampleSet(const SegmentedSampleSet<V>& set)` is currently a stub
  - currently throws even when marker information is present
  - requires defining behavior for markers not covered by any segment

## Cleaning features

- [ ] Implement cleaning for raw CN files in `src/cna_clean.hpp`
  - `Clean::clean(segmented<false>)` currently throws

## Sample / container behavior

- [ ] Implement `Sample::remove(T)` in `lib/Sample.hpp`
  - currently throws `"Sample::remove(T) has yet been implemented."`

## Default value interference

- [ ] Determine sane default values for stateDiff and refState for cna clean and
      filter based on input type (current defaults are for LRR input)

## Chromosome number

- [ ] Support different number of chromosomes (i.e. for non-humans)

