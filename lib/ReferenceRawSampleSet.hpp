#ifndef cna_ReferenceRawSampleSet_h
#define cna_ReferenceRawSampleSet_h

namespace cna {

// same as RawSampleSet
template <typename V>
class ReferenceRawSampleSet : public RawSampleSet<V>
{
public:
	typedef RawSampleSet<V> Base;
};

} // namespace cna


#endif
