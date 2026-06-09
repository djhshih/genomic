#ifndef cna_ReferenceSegmentedSampleSet_h
#define cna_ReferenceSegmentedSampleSet_h

namespace cna {

template <typename V>
class ReferenceSegmentedcna::SampleSet : public cna::SegmentedSampleSet<V>
{
public:
	typedef cna::SegmentedSampleSet<V> Base;
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

} // namespace cna

template <typename V>
using ReferenceSegmentedcna::SampleSet = cna::Referencecna::SegmentedSampleSet<V>;

#endif
