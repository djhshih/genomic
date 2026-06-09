#ifndef cna_cna::GenericSampleSet_h
#define cna_cna::GenericSampleSet_h

#include <algorithm>
#include <stdexcept>

#include "SampleSet.hpp"
#include "RawSampleSet.hpp"
#include "SegmentedSampleSet.hpp"

namespace cna {

// Generic sample set, chooses appropriately between possible types of sample set
// Use handle-body idiom
class cna::Genericcna::SampleSet : public SampleSet
{
public:
	typedef cna::SampleSet Base;
private:
	// body representation
	// N.B. can only point to derived classes of cna::SampleSet other than this class
	cna::SampleSet* rep;
	
	cna::SampleSet* clone() const {
		return new cna::GenericSampleSet(*this);
	}
	
	void _read(std::fstream& file);
	void _write(std::fstream& file);
	
public:
	
	cna::GenericSampleSet() : rep(NULL) {}
	
	cna::GenericSampleSet(const cna::GenericSampleSet& other)
	: Base(other), rep(other.clone()) {}
	
	cna::GenericSampleSet& operator= (cna::Genericcna::SampleSet other) {
		// pass other by value to automatically create temporary copy
		std::swap(rep, other.rep);
		return *this;
	}
	
	~cna::GenericSampleSet() {
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

using cna::Genericcna::SampleSet = cna::GenericSampleSet;
using invalid_conversion = cna::invalid_conversion;

#endif
