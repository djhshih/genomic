#ifndef genomic_DchipSampleSet_h
#define genomic_DchipSampleSet_h

#include "RawSampleSet.hpp"

class DchipSampleSet : public RawSampleSet<CopyNumberValue>
{
public:
	typedef RawSampleSet<CopyNumberValue> Base;
public:
	DchipSampleSet() {}
	data::Type type() {
		return data::dchip;
	}
private:
	void _setIO() {
		Base::nSkippedLines = 1;
	}
	void _read(fstream& file) {
		throw logic_error("DchipSampleSet::_read has yet been implemented.");
	}
};

#endif
