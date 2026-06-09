#ifndef cna_ReferenceRawSampleSet_h
#define cna_ReferenceRawSampleSet_h

namespace cna {

// same as RawSampleSet
template <typename V>
class ReferenceRawcna::SampleSet : public cna::RawSampleSet<V>
{
public:
	typedef cna::RawSampleSet<V> Base;
};

} // namespace cna

template <typename V>
using ReferenceRawcna::SampleSet = cna::Referencecna::RawSampleSet<V>;

#endif
