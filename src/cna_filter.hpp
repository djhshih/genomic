#ifndef cna_filter_h
#define cna_filter_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "cna_common.hpp"


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
			cna::dice_overlapper checker(threshold);
			set.template filter<cna::dice_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "query") {
			cna::query_overlapper checker(threshold);
			set.template filter<cna::query_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "reference") {
			cna::reference_overlapper checker(threshold);
			set.template filter<cna::reference_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "min") {
			cna::min_overlapper checker(threshold);
			set.template filter<cna::min_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else if (score == "max") {
			cna::max_overlapper checker(threshold);
			set.template filter<cna::max_overlapper>(ref, checker, inverse, merge, aberrant, optimize);
		} else {
			throw std::runtime_error("Invalid overlap score method specified.");
		}
		
		set.write(outputFileName);
	}
	
	void run();
	
private:
	
	std::string inputFileName, referenceFileName, outputFileName;
	cna::data::Type inputType, referenceType;
	float threshold;
	bool merge, aberrant;
	float stateDiff, refState;
	bool optimize;
	bool inverse;
	std::string score;
	
	void getOptions();
	
};

#endif
