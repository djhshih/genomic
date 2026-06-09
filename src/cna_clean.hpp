#ifndef cna_clean_h
#define cna_clean_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "cna_common.hpp"


class Clean : public Command {
	
public:
	
	Clean();
	
	template <bool> struct segmented {};
	
	template <typename SampleSetType>
	void clean(segmented<false>) {
		throw std::logic_error("Cleaning functions for raw CN files have yet been implemented.");
	}
	
	template <typename SampleSetType>
	void clean(segmented<true>) {
		SampleSetType set;
		set.read(inputFileName);
		
		set.set(CNACriteria(refState, stateDiff));
		
		if (count > 0) set.filter(cna::spurious_segment_filter<typename SampleSetType::Value>(count), inverse, merge);
		if (length > 0) set.filter(cna::small_segment_filter<typename SampleSetType::Value>(length), inverse, merge);
		if (balanced) set.filter(cna::balanced_segment_filter<typename SampleSetType::Value>(refState, stateDiff), inverse, false);
		
		set.write(outputFileName);
	}
	
	void run();
	
private:
	
	std::string inputFileName, outputFileName;
	cna::data::Type inputType;
	float diceThreshold;
	bool inverse, merge, balanced;
	position count, length;
	float stateDiff, refState;
	
	void getOptions();
	
};

#endif
