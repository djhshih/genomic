#ifndef genomic_PicnicSampleSet_h
#define genomic_PicnicSampleSet_h

class PicnicSampleSet : public SplitRawSampleSet<AlleleSpecificIntegerCopyNumberValue>
{
public:
	typedef SplitRawSampleSet<AlleleSpecificIntegerCopyNumberValue> Base;
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
		Base::Base::delim = ',';
		Base::Base::headerLine = 0;
	}
};

#endif
