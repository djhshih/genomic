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
	
	Filter() 
	: Command("filter file based on reference") {
		
		// Delcare options
		opts.add_options()
			("help", "print help message")
			("input,i", po::value<string>(), "sample set file")
			("reference,r", po::value<string>(), "reference set file (same format as sample set)")
			("format,f", po::value<string>(), "input file format [default: determined from file extension]")
			("reference_format,g", po::value<string>(), "reference file format [default: reference format corresponding to input]")
			("output,o", po::value<string>(), "output file")
			("threshold,t", po::value<float>(), "overlap threshold (only used for segmentation files) [default: 0.5]")
			("score,s", po::value<string>(), "overlap score method, only used for segmentation files) [options: dice (default), query, reference, min, max]")
			("inverse,v", po::value<bool>(), "select instead of filter overlapping segments")
			("merge,m", po::value<bool>(), "merge filtered segments with upstream/downstream segments?")
			("aberrant,a", po::value<bool>(), "filter only aberrant segments?")
			("state_diff", po::value<rvalue>(), "threshold for difference from reference state")
			("ref_state", po::value<rvalue>(), "reference state")
			("optimize,O", po::value<bool>(), "optimize algorithm speed, assuming contiguity of reference segments")
			;
		popts.add("input", 1).add("reference", 1).add("output", 1);
		
	}
	
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
			throw runtime_error("Invalid overlap score method specified.");
		}
		
		set.write(outputFileName);
	}
	
	void run() {
		
		if (vm.count("help")) {
			cout << "usage:  " << progname << " filter [options] <sample set file> <reference set file> <output file>" << endl;
			cout << opts << endl;
			return;
		}
		
		getOptions();
		
		switch (inputType) {
			
			case data::segmented:
				
				switch(referenceType) {
					case data::segmented_ref:
						filter< SegmentedSampleSet<rvalue>, ReferenceSegmentedSampleSet<rvalue> >(segmented<true>());
						break;
					case data::segmented:
						filter< SegmentedSampleSet<rvalue>, SegmentedSampleSet<rvalue> >(segmented<true>());
						break;
					default:
						throw invalid_argument("Unsupported reference file format.");
				}
				
				break;
				
			case data::segmented_ascn:
				
				switch(referenceType) {
					case data::segmented_ref:
						filter< SegmentedSampleSet<alleles_cn>, ReferenceSegmentedSampleSet<alleles_cn> >(segmented<true>());
						break;
					case data::segmented_ascn:
						filter< SegmentedSampleSet<alleles_cn>, SegmentedSampleSet<alleles_cn> >(segmented<true>());
						break;
					default:
						throw invalid_argument("Unsupported reference file format.");
				}
				
				break;
				
			case data::raw:
				
				switch(referenceType) {
					case data::raw_ref:
						filter< RawSampleSet<rvalue>, ReferenceRawSampleSet<rvalue> >(segmented<false>());
						break;
					case data::raw:
						filter< RawSampleSet<rvalue>, RawSampleSet<rvalue> >(segmented<false>());
						break;
					default:
						throw invalid_argument("Unsupported reference file format.");
				}
				
				break;
				
			case data::raw_ascn:
				
				switch(referenceType) {
					case data::raw_ref:
						filter< RawSampleSet<alleles_cn>, ReferenceRawSampleSet<alleles_cn> >(segmented<false>());
						break;
					case data::raw_ascn:
						filter< RawSampleSet<alleles_cn>, RawSampleSet<alleles_cn> >(segmented<false>());
						break;
					default:
						throw invalid_argument("Unsupported reference file format.");
				}
				
				break;
				
		}
		
		
	}
	
private:
	
	string inputFileName, referenceFileName, outputFileName;
	data::Type inputType, referenceType;
	float threshold;
	bool merge, aberrant;
	float stateDiff, refState;
	bool optimize;
	bool inverse;
	string score;
	
	void getOptions() {
		
		if (vm.count("input") && vm.count("reference")) {
			inputFileName = vm["input"].as<string>();
			referenceFileName = vm["reference"].as<string>();
		} else {
			throw invalid_argument("Both input and reference files must be specified.");
		}
		
		if (vm.count("format")) {
			inputType = mapping::extension[ vm["format"].as<string>() ];
		} else {
			inputType = mapping::extension[ name::fileext(inputFileName) ];
		}
		if (inputType == data::invalid) {
			throw invalid_argument("Invalid input format type.");
		}
		
		if (vm.count("reference_format")) {
			referenceType = mapping::extension[ vm["reference_format"].as<string>() ];
		} else {
			if (inputType == data::raw || inputType == data::raw_ascn || inputType == data::raw_lrrbaf) {
				referenceType = data::raw_ref;
			} else {
				referenceType = data::segmented_ref;
			}
		}
		if (referenceType == data::invalid) {
			throw invalid_argument("Invalid reference format type.");
		}
		
		if (vm.count("output")) {
			outputFileName = vm["output"].as<string>();
		} else {
			outputFileName = name::filestem(inputFileName) + ".filtered." + name::fileext(inputFileName);
		}
		
		if (vm.count("threshold")) {
			threshold = vm["threshold"].as<float>();
		} else {
			threshold = 0.5;
		}
		
		if (vm.count("score")) {
			score = vm["score"].as<string>();
		} else {
			score = "dice";
		}
		
		if (vm.count("merge")) {
			merge = vm["merge"].as<bool>();
		} else {
			merge = false;
		}
		
		if (vm.count("aberrant")) {
			aberrant = vm["aberrant"].as<bool>();
		} else {
			aberrant = false;
		}
		
		//TODO automatically determine stateDiff and refState
		// current settings are suitable for LRR data
		
		if (vm.count("state_diff")) {
			stateDiff = vm["state_diff"].as<float>();
		} else {
			stateDiff = 0.1;
		}
		
		if (vm.count("ref_state")) {
			refState = vm["ref_state"].as<float>();
		} else {
			refState = 0;
		}
		
		if (vm.count("optimize")) {
			optimize = vm["optimize"].as<bool>();
		} else {
			optimize = true;
		}
		
		if (vm.count("inverse")) {
			inverse = vm["inverse"].as<bool>();
		} else {
			inverse = false;
		}
	}
	
};

#endif
