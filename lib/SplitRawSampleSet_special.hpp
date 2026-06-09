// define SPECIALIZATION_TYPE
// include in SplitRawSampleSet.h multiple times with different SPECIALIZATION_TYPE

template <> inline
void cna::SplitRawSampleSet<SPECIALIZATION_TYPE>::readSampleValue(std::istringstream& stream, typename Base::RawSample* sample, size_t chromIndex, const char delim) {
	typename Base::Value value;
	std::string s;
	getline(stream, s, delim);
	value.a = atof(s.c_str());
	getline(stream, s, delim);
	value.b = atof(s.c_str());
	sample->addToChromosome(chromIndex, value);
}
