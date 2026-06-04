#ifndef genomic_PicnicSampleSet_h
#define genomic_PicnicSampleSet_h

#include "typedefs.h"

template <typename V> class SplitRawSampleSet;

class PicnicSampleSet : public SplitRawSampleSet<alleles_cn>
{
public:
	typedef SplitRawSampleSet<alleles_cn> Base;
public:
	// set dataColumn to 5
	PicnicSampleSet() : Base(5) {
		_setIO();
	}
	data::Type type() {
		return data::picnic;
	}
private:
	void _setIO() {
		Base::Base::io.delim = ',';
		Base::Base::io.headerLine = 0;
	}
};

#endif
