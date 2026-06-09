#ifndef cna_GenericSampleSet_h
#define cna_GenericSampleSet_h

#include <algorithm>
#include <stdexcept>

#include "SampleSet.hpp"
#include "RawSampleSet.hpp"
#include "SegmentedSampleSet.hpp"

namespace cna {

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
	
	SampleSet* clone() const {
		return new GenericSampleSet(*this);
	}
	
	void _read(std::fstream& file);
	void _write(std::fstream& file);
	
public:
	
	GenericSampleSet() : rep(NULL) {}
	
	GenericSampleSet(const GenericSampleSet& other)
	: Base(other), rep(other.clone()) {}
	
	GenericSampleSet& operator= (GenericSampleSet other) {
		// pass other by value to automatically create temporary copy
		std::swap(rep, other.rep);
		return *this;
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
	
	cna::data::Type type() {
		return cna::data::generic;
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
	explicit invalid_conversion(const std::string& what_arg)
	: std::logic_error("Invalid conversion. " + what_arg) {}
};

} // namespace cna


#endif
