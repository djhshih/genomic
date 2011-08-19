// define SPECIALIZATION_TYPE
// include in SplitRawSampleSet.h multiple times with different SPECIALIZATION_TYPE

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
