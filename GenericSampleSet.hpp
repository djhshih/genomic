#ifndef genomic_GenericSampleSet_h
#define genomic_GenericSampleSet

#include <stdexcept>

#include "SampleSet.hpp"
#include "RawSampleSet.hpp"
#include "SegmentedSampleSet.hpp"

//template <typename V> class SegmentedSampleSet;
//template <typename V> class RawSampleSet;

// Generic sample set, chooses appropriately between possible types of sample set
// Use handle-body idiom
class GenericSampleSet : public SampleSet
{
public:
	typedef SampleSet Base;
private:
	// body representation
	// N.B. can only point to derived classes of SampleSet other than this class
	SampleSet* rep;
	
	SampleSet* clone() {
		return new GenericSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);
public:
	GenericSampleSet() : rep(NULL) {}
	GenericSampleSet(const GenericSampleSet& gset) {
		rep = gset.rep->clone();
	}
	~GenericSampleSet() {
		clear();
	}
	void clear() {
		delete rep;
		rep = NULL;
	}
	void sort() {
		if (rep != NULL) {
			rep->sort();
		}
	}
	data::Type type() {
		return data::generic;
	}
	size_t size() {
		if (rep != NULL) {
			return rep->size();
		}
		return 0;
	}
};

class invalid_conversion : public std::logic_error
{
public:
	explicit invalid_conversion(const string& what_arg)
	: std::logic_error("Invalid conversion. " + what_arg) {}
};

#endif
