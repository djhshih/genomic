# Architecture Overview

## Purpose
`genomic` is a C++ command-line toolkit for reading, transforming, filtering, cleaning, and sorting genomic copy-number style data. The codebase is organized around a small CLI front end and a reusable in-memory data model for markers, samples, chromosomes, raw measurements, and segmented intervals.

## High-level structure

### 1. Application entry point
- **`genomic.cpp`**
  - Creates the available subcommands: `convert`, `filter`, `clean`, and `sort`.
  - Dispatches based on the first CLI argument.
  - Prints usage/version information.

### 2. Command layer
- **`genomic_common.hpp/cpp`**
  - Defines the `Command` base class used by all CLI commands.
  - Wraps Boost Program Options parsing.
  - Exposes global `progname` and stream formatting for commands.
- **`genomic_convert.hpp`**
  - Implements format conversion workflows.
- **`genomic_filter.hpp`**
  - Implements filtering against reference datasets.
- **`genomic_clean.hpp`**
  - Implements rule-based cleanup of segmented data.
- **`genomic_sort.hpp`**
  - Implements marker-file sorting.

This layer is thin: it mostly parses options, chooses the right dataset types, and invokes domain operations.

## Core domain model

### 3. Global types, constants, and helpers
- **`typedefs.h`**
  - Defines core scalar aliases used across the project.
- **`global.hpp/cpp`**
  - Global chromosome constants (`nAutosomes`, `nChromosomes`).
  - Utility templates like `absdiff`, `eq`, `neq`.
  - Name/path helpers in namespace `name`.
  - Lookup singletons in namespace `mapping`:
    - chromosome string ↔ numeric id mapping
    - file extension ↔ `data::Type` mapping
- **`Properties.hpp`**
  - `IOProperties` for parsing/writing behavior.
  - `CNACriteria` for copy-number state thresholds.

These files act as shared infrastructure for the rest of the system.

### 4. Generic biological/container hierarchy
The domain model is layered from low-level genomic coordinates upward:

- **Markers** describe probe/marker positions on chromosomes.
- **Chromosomes** hold ordered sequences of values or segments.
- **Samples** group chromosome data for one biological sample.
- **Sample sets** group many samples and own the file I/O behavior.

#### `Marker` subsystem
- **`Marker.hpp/cpp`**
  - `marker::Marker` stores marker name, chromosome, genomic position, and filter flag.
  - `marker::Set` stores per-chromosome collections of markers.
  - Supports reading, writing, sorting, filtering, distribution of unsorted markers, and cleanup.
  - `marker::Manager` is a global registry/reference-count manager for shared marker sets.

This is a key architectural choice: multiple sample sets can share the same marker set by platform name instead of duplicating marker metadata.

#### `Sample` abstraction
- **`Sample.hpp`**
  - Templated `Sample<Chromosome>`.
  - Holds one chromosome container per chromosome.
  - Provides indexed and name-based chromosome access.
  - Used by both raw and segmented representations.

#### `SampleSet` abstraction
- **`SampleSet.hpp/cpp`**
  - Abstract base class for all dataset containers.
  - Owns common read/write orchestration and file handles.
  - Stores shared state:
    - `IOProperties`
    - `CNACriteria`
    - associated `marker::Set*`
    - source/output filename
  - Defines virtual hooks:
    - `_read`
    - `_write`
    - `clear`
    - `sort`
    - `type`
    - `size`
    - `clone`

Important design point: `SampleSet` handles the generic workflow of opening files and coordinating marker setup, while concrete subclasses implement the actual serialization logic.

## Concrete dataset representations

### 5. Raw data representation
- **`RawSampleSet.hpp`**
  - Template `RawSampleSet<V>`.
  - Stores per-marker values for each sample.
  - Internally uses:
    - `Sample<LinearChromosome<Value>>`
    - vectors of sample pointers
    - sample-name lookup map
  - Responsibilities:
    - parse matrix-style raw CN files
    - write raw CN files
    - sort raw values consistently with markers
    - filter by marker membership using a reference marker set

The raw representation is marker-aligned: values are stored in the same order as the associated marker set.

### 6. Segmented data representation
- **`SegmentedSampleSet.hpp`**
  - Template `SegmentedSampleSet<V>`.
  - Stores interval segments rather than marker-level values.
  - Internally uses `Sample<LinearChromosome<Segment<Value>>>`.
  - Responsibilities:
    - parse segment tables
    - write segment tables
    - sort samples and segments
    - convert from raw data by collapsing adjacent equal-valued markers into segments
    - perform several filtering/cleanup operations

This file also contains much of the domain logic for segment processing:
- filter functor interface (`filter_operator`)
- concrete filters such as:
  - `spurious_segment_filter`
  - `small_segment_filter`
  - `balanced_segment_filter`
  - `reference_segment_filter`
- overlap scoring strategies:
  - dice
  - query
  - reference
  - min/max
- merging and removal logic for flagged segments

Architecturally, this is the richest part of the codebase: the segmented representation is where most analytical behavior lives.

### 7. Dynamic wrapper for format-driven behavior
- **`GenericSampleSet.hpp/cpp`**
  - Implements a handle/body style wrapper over `SampleSet*`.
  - Determines the concrete representation from file extension at read/write time.
  - Delegates to `RawSampleSet` or `SegmentedSampleSet` variants.
  - Performs limited runtime conversion between raw and segmented forms.

This gives the CLI a type-erased entry point when the exact file type is not known at compile time.

### 8. Format-specific specializations
- **`ReferenceRawSampleSet.hpp`**
- **`ReferenceSegmentedSampleSet.hpp`**
- **`PicnicSampleSet.hpp`**
- **`CnagSampleSet.hpp`**
- **`DchipSampleSet.hpp`**
- **`SplitRawSampleSet.hpp`**
- **`RawSampleSet_special.hpp`**
- **`SegmentedSampleSet_special.hpp`**

These files extend the core raw/segmented model for particular external formats or value types, especially allele-specific forms. The main architecture remains the same: they fit under the `SampleSet` family and reuse shared parsing, storage, and conversion patterns.

## Data flow

## 9. Typical execution flow

### Convert
1. CLI dispatch selects `Convert`.
2. Command determines input/output `data::Type` from extension or flags.
3. Appropriate concrete sample set is instantiated.
4. Input is read into in-memory structures.
5. Optional conversion between raw and segmented forms occurs.
6. Output is written in the target format.

### Filter
1. CLI dispatch selects `Filter`.
2. Sample and reference datasets are loaded as matching concrete types.
3. For raw data, marker membership filtering is applied.
4. For segmented data, overlap-based filters compare segments to the reference.
5. Flagged segments/markers are removed or optionally merged.
6. Result is written out.

### Clean
1. Segmented dataset is loaded.
2. One or more filter functors are applied.
3. Flagged segments are removed or merged.
4. Cleaned result is written out.

### Sort
1. Marker file is read.
2. Markers are sorted by chromosome and genomic position.
3. Sorted marker file is written.

## Architectural patterns used

### 10. Main patterns
- **Template-based polymorphism**
  - Used heavily for value types (`rvalue`, allele-specific values, segments).
- **Runtime polymorphism via base class**
  - `SampleSet` and `Command` provide virtual interfaces.
- **Handle/body type erasure**
  - `GenericSampleSet` hides concrete dataset type behind `SampleSet*`.
- **Strategy/functor pattern**
  - Segment filtering and overlap scoring are implemented as interchangeable functors/classes.
- **Shared resource manager**
  - `marker::Manager` centralizes marker-set lifetime and reuse.

## Module boundaries

### 11. Main dependencies between modules
- `genomic.cpp` depends on the command classes.
- Command classes depend on:
  - `genomic_common`
  - `global`
  - `SampleSets.hpp`
- `SampleSets.hpp` aggregates all sample-set variants.
- `SampleSet` depends on `Sample`, `Marker`, and `Properties`.
- `RawSampleSet` and `SegmentedSampleSet` derive from `SampleSet` and contain most format/domain logic.
- `Marker` depends on global chromosome mapping and IO settings.

A simplified dependency graph is:

```text
CLI main
  -> Command subclasses
    -> SampleSet family
      -> Sample / Chromosome / Segment / Marker
        -> global typedefs + mappings + properties
```

## Build and test structure

### 12. Build
- **`CMakeLists.txt`** builds the main executable from a small set of translation units:
  - `genomic.cpp`
  - `genomic_common.cpp`
  - `global.cpp`
  - `SampleSet.cpp`
  - `GenericSampleSet.cpp`
  - `Marker.cpp`
- Much of the functional code is header-only through templates.

### 13. Tests
- **`tests/genomic_test.cpp`**
  - Exercises I/O, conversions, copy behavior, filtering, and mapping utilities.
  - Confirms the architecture’s main contract: datasets can be read, transformed, and written reproducibly.

## Summary
The codebase is centered on a reusable in-memory model for genomic sample data, with a CLI layer on top. The dominant architectural idea is a combination of:
- abstract base classes for shared workflows,
- templates for value/type specialization,
- and file-extension-driven dispatch for format selection.

In practice, the most important subsystems are:
1. **CLI commands** for user-facing workflows
2. **`SampleSet` hierarchy** for data loading/writing/conversion
3. **marker management** for shared coordinate metadata
4. **segmented-data filtering logic** for the core analytical operations

The graphics code exists as a separate experimental branch of the architecture, but the primary product is the CLI data-processing pipeline.
