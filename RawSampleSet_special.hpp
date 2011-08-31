// define SPECIALIZATION_TYPE
// include in RawSampleSet.hpp multiple times with different SPECIALIZATION_TYPE

template <> inline
data::Type RawSampleSet<SPECIALIZATION_TYPE>::type() {
	return data::raw_ascn;
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::readSampleNames(istringstream& stream) {
	while (!stream.eof()) {
		string sampleName1, sampleName2;
		stream >> sampleName1 >> sampleName2;
		create(name::common(sampleName1, sampleName2));
	}
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::readSampleValues(istringstream& stream, size_t sampleStart, const string& chromName) {
	size_t i = sampleStart;
	Value value;
	while (!stream.eof()) {
		stream >> value.a >> value.b;
		// create point at specified chromosome
		samples[++i]->addToChromosome(chromName, value);
	}
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::writeSampleNames(fstream& file, const char delim) {
	// print sample names
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name << ".A" << delim << (**it).name << ".B";
	}
	file << endl;
}

template <> inline
void RawSampleSet<SPECIALIZATION_TYPE>::writeSampleValues(fstream& file, size_t chr, size_t markerIndex, const char delim) {
	// iterate through samples to print values, selected the specified chromosome and marker
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		const Value& value = (**it)[chr][markerIndex];
		file << delim << value.a << delim << value.b;
	}
	file << endl;
}
