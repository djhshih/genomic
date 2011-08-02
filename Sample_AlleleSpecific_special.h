// define SPECIALIZATION_TYPE
// include in Sample.h

template <> inline
data::Type RawSampleSet<SPECIALIZATION_TYPE>::type() {
	return data::raw_ascn;
}

template <> inline
data::Type SegmentedSampleSet<SPECIALIZATION_TYPE>::type() {
	return data::segmented_ascn;
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
void SegmentedSampleSet<SPECIALIZATION_TYPE>::readSegment(fstream& file, Segment<SPECIALIZATION_TYPE>& seg) {
	file >> seg.start >> seg.end;
	if (!positionsOnly) {
		file >> seg.count >> seg.value.a >> seg.value.b;
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
		const Value& value = (**it)[chr]->at(markerIndex);
		file << delim << value.a << delim << value.b;
	}
	file << endl;
}

template <> inline
void SegmentedSampleSet<SPECIALIZATION_TYPE>::_write(fstream& file)
{
	const char delim = Base::delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "stateA" << delim << "stateB" << endl;
	
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->count << delim << segIt->value.a << delim << segIt->value.b << endl;
			}
			++chr;
		}
	}
}

template <> inline
void SplitRawSampleSet<SPECIALIZATION_TYPE>::readSampleValue(istringstream& stream, typename Base::RawSample* sample, size_t chromIndex, const char delim) {
	typename Base::Value value;
	string s;
	getline(stream, s, delim);
	value.a = atof(s.c_str());
	getline(stream, s, delim);
	value.b = atof(s.c_str());
	sample->addToChromosome(chromIndex, value);
}
