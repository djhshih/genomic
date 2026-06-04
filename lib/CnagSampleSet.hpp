#ifndef genomic_CnagSampleSet_h
#define genomic_CnagSampleSet_h

#include "RawSampleSet.hpp"

class CnagSampleSet : public RawSampleSet<cnvalue>
{
public:
	typedef RawSampleSet<cnvalue> Base;
public:
	CnagSampleSet() {}
	data::Type type() {
		return data::cnag;
	}
private:
	//void _setIO() {}
	void _read(fstream& file) {
		throw logic_error("CnagSampleSet::_read has yet been implemented.");
	}
};

#endif
