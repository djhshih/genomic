#ifndef genomic_ReferenceRawSampleSet_h
#define genomic_ReferenceRawSampleSet_h

template <typename V> class RawSampleSet;

// same as RawSampleSet
template <typename V>
class ReferenceRawSampleSet : public RawSampleSet<V>
{
public:
	typedef RawSampleSet<V> Base;
};

#endif
