#ifndef genomic_filter_h
#define genomic_filter_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "genomic_common.hpp"


class Filter : public Command {
	
public:
	
	Filter();
	
	template <bool> struct segmented {};
	
	template <typename SampleSetType, typename ReferenceSetType>
	void filter(segmented<false>) {
		SampleSetType set;
		set.read(inputFileName);
		ReferenceSetType ref;
		ref.read(referenceFileName);
		set.filter(ref);
		set.write(outputFileName);
	}
	
	template <typename SampleSetType, typename ReferenceSetType>
	void filter(segmented<true>) {
		SampleSetType set;
		set.read(inputFileName);
		
		ReferenceSetType ref;
		ref.read(referenceFileName);
		
		set.set(CNACriteria(refState, stateDiff));
		
		if (score == "dice") {
			dice_overlapper checker(threshold);
			set.template filter<dice_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "query") {
			query_overlapper checker(threshold);
			set.template filter<query_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "reference") {
			reference_overlapper checker(threshold);
			set.template filter<reference_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "min") {
			min_overlapper checker(threshold);
			set.template filter<min_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "max") {
			max_overlapper checker(threshold);
			set.template filter<max_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else {
			throw std::runtime_error("Invalid overlap score method specified.");
		}
		
		set.write(outputFileName);
	}
	
	void run();
	
private:
	
	std::string inputFileName, referenceFileName, outputFileName;
	data::Type inputType, referenceType;
	float threshold;
	bool merge, aberrant;
	float stateDiff, refState;
	bool optimize;
	bool inverse;
	std::string score;
	
	void getOptions();
	
};

#endif
