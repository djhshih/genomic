#ifndef genomic_clean_h
#define genomic_clean_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "genomic_common.hpp"


class Clean : public Command {
	
public:
	
	Clean() 
	: Command("clean file based on filtering criteria") {
		
		// Delcare options
		opts.add_options()
			("help", "print help message")
			("input,i", po::value<string>(), "sample set file")
			("format,f", po::value<string>(), "input file format [default: determined from file extension]")
			("output,o", po::value<string>(), "output file")
			("merge,m", po::value<bool>(), "merge filtered segments with upstream/downstream segments?")
			("count", po::value<position>(), "threshold for number of markers in segment")
			("length", po::value<position>(), "threshold for segment length")
			("balanced", po::value<bool>(), "remove balanced segments?")
			("state_diff", po::value<rvalue>(), "threshold for difference from reference state")
			("ref_state", po::value<rvalue>(), "reference state")
			;
		popts.add("input", 1).add("reference", 1).add("output", 1);
		
	}
	
	template <bool> struct segmented {};
	
	template <typename SampleSetType>
	void clean(segmented<false>) {
		throw logic_error("Cleaning functions for raw CN files have yet been implemented.");
	}
	
	template <typename SampleSetType>
	void clean(segmented<true>) {
		SampleSetType set;
		set.read(inputFileName);
		
		set.set(CNACriteria(refState, stateDiff));
		
		if (count > 0) set.filter(spurious_segment_filter<typename SampleSetType::Value>(count), merge);
		if (length > 0) set.filter(small_segment_filter<typename SampleSetType::Value>(length), merge);
		if (balanced) set.filter(balanced_segment_filter<typename SampleSetType::Value>(refState, stateDiff), false);
		
		set.write(outputFileName);
	}
	
	void run() {
		
		if (vm.count("help")) {
			cout << "usage:  " << progname << " clean [options] <sample set file> <reference set file> <output file>" << endl;
			cout << opts << endl;
			return;
		}
		
		getOptions();
		
		switch (inputType) {
			
			case data::segmented:
				
				clean< SegmentedSampleSet<rvalue> >(segmented<true>());
				break;
				
			case data::segmented_ascn:
				
				clean< SegmentedSampleSet<alleles_cn> >(segmented<true>());
				break;
				
			case data::raw:
				
				clean< RawSampleSet<rvalue> >(segmented<false>());
				break;
				
			case data::raw_ascn:
				
				clean< RawSampleSet<alleles_cn> >(segmented<false>());
				break;
				
		}
		
	}
	
private:
	
	string inputFileName, outputFileName;
	data::Type inputType;
	float diceThreshold;
	bool merge, balanced;
	position count, length;
	float stateDiff, refState;
	
	void getOptions() {
		
		if (vm.count("input")) {
			inputFileName = vm["input"].as<string>();
		} else {
			throw invalid_argument("Input file not specified.");
		}
		
		if (vm.count("format")) {
			inputType = mapping::extension[ vm["format"].as<string>() ];
		} else {
			inputType = mapping::extension[ name::fileext(inputFileName) ];
		}
		if (inputType == data::invalid) {
			throw invalid_argument("Invalid input format type.");
		}
		
		if (vm.count("output")) {
			outputFileName = vm["output"].as<string>();
		} else {
			outputFileName = name::filestem(inputFileName) + ".filtered." + name::fileext(inputFileName);
		}
		
		if (vm.count("merge")) {
			merge = vm["merge"].as<bool>();
		} else {
			merge = false;
		}
		
		if (vm.count("count")) {
			count = vm["count"].as<position>();
		} else {
			count = 0;
		}
		
		if (vm.count("length")) {
			length = vm["length"].as<position>();
		} else {
			length = 0;
		}
		
		//TODO automatically determine stateDiff and refState
		// current settings are suitable for LRR data
		
		if (vm.count("balanced")) {
			balanced = vm["balanced"].as<bool>();
		} else {
			balanced = false;
		}
		
		if (vm.count("state_diff")) {
			stateDiff = vm["state_diff"].as<float>();
		} else {
			stateDiff = 0.2;
		}
		
		if (vm.count("ref_state")) {
			refState = vm["ref_state"].as<float>();
		} else {
			refState = 0;
		}
	}
	
};

#endif
