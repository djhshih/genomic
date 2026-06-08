#ifndef cna_ReferenceSegmentedSampleSet_h
#define cna_ReferenceSegmentedSampleSet_h

template <typename V> class SegmentedSampleSet;

template <typename V>
class ReferenceSegmentedSampleSet : public SegmentedSampleSet<V>
{
public:
	typedef SegmentedSampleSet<V> Base;
public:
	ReferenceSegmentedSampleSet() {
		_setIO();
	}
private:
	void _setIO() {
		Base::mergeSamples = true;
		Base::positionsOnly = true;
	}
};

#endif
