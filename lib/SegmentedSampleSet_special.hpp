// define SPECIALIZATION_TYPE
// include in SegmentedSampleSet.h multiple times with different SPECIALIZATION_TYPE

template <> inline
cna::data::Type SegmentedSampleSet<SPECIALIZATION_TYPE>::type() {
	return cna::data::segmented_ascn;
}

template <> inline
void SegmentedSampleSet<SPECIALIZATION_TYPE>::readSegment(FieldScanner& fields, Segment<SPECIALIZATION_TYPE>& seg) {
	std::string_view field;
	if (!fields.next(field) || !parseNumber(field, seg.start)) return;
	if (!fields.next(field) || !parseNumber(field, seg.end)) return;
	if (!positionsOnly) {
		if (!fields.next(field) || !parseNumber(field, seg.count)) return;
		if (!fields.next(field) || !parseNumber(field, seg.value.a)) return;
		if (!fields.next(field) || !parseNumber(field, seg.value.b)) return;
	}
}

template <> inline
void SegmentedSampleSet<SPECIALIZATION_TYPE>::_write(std::fstream& file)
{
	const char delim = Base::io.delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "stateA" << delim << "stateB" << std::endl;
	
	typename Samples::iterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		typename Chromosomes::iterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			typename Segments::iterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << segIt->chromosome << delim << segIt->start << delim << segIt->end << delim << segIt->count << delim << segIt->value.a << delim << segIt->value.b << std::endl;
			}
		}
	}
}
