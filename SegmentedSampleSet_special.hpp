// define SPECIALIZATION_TYPE
// include in SegmentedSampleSet.h multiple times with different SPECIALIZATION_TYPE

template <> inline
data::Type SegmentedSampleSet<SPECIALIZATION_TYPE>::type() {
	return data::segmented_ascn;
}

template <> inline
void SegmentedSampleSet<SPECIALIZATION_TYPE>::readSegment(istringstream& stream, Segment<SPECIALIZATION_TYPE>& seg) {
	stream >> seg.start >> seg.end;
	if (!positionsOnly) {
		stream >> seg.count >> seg.value.a >> seg.value.b;
	}
}

template <> inline
void SegmentedSampleSet<SPECIALIZATION_TYPE>::_write(fstream& file)
{
	const char delim = Base::io.delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "stateA" << delim << "stateB" << endl;
	
	typename Samples::iterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		typename Chromosomes::iterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			typename Segments::iterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << segIt->chromosome << delim << segIt->start << delim << segIt->end << delim << segIt->count << delim << segIt->value.a << delim << segIt->value.b << endl;
			}
		}
	}
}
