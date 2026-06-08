// define SPECIALIZATION_TYPE
// include in RawSampleSet.hpp multiple times with different SPECIALIZATION_TYPE

template <> inline
data::Type RawSampleSet<SPECIALIZATION_TYPE>::type() {
	return data::raw_ascn;
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::readSampleNames(FieldScanner& fields) {
	std::string_view sampleName1, sampleName2;
	while (fields.next(sampleName1) && fields.next(sampleName2)) {
		create(name::common(std::string(sampleName1), std::string(sampleName2)));
	}
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::readSampleValues(FieldScanner& fields, size_t sampleStart, const std::string& chromName) {
	size_t i = sampleStart;
	Value value;
	std::string_view a, b;
	while (fields.next(a) && fields.next(b)) {
		if (!parseNumber(a, value.a) || !parseNumber(b, value.b)) {
			continue;
		}
		// create point at specified chromosome
		samples[++i]->addToChromosome(chromName, value);
	}
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::writeSampleNames(std::fstream& file, const char delim) {
	// print sample names
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name << ".A" << delim << (**it).name << ".B";
	}
	file << std::endl;
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::writeSampleValues(std::fstream& file, size_t chr, size_t markerIndex, const char delim) {
	// iterate through samples to print values, selected the specified chromosome and marker
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		const Value& value = (**it)[chr][markerIndex];
		file << delim << value.a << delim << value.b;
	}
	file << std::endl;
}
