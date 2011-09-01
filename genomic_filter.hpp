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
			("dice,d", po::value<float>(), "Dice coefficient threshold (only used for segmentation files) [default: 0.8]")
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
		set.filter(ref, diceThreshold);
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
	float diceThreshold;
	
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
		
		if (vm.count("dice")) {
			diceThreshold = vm["dice"].as<float>();
		} else {
			diceThreshold = 0.8;
		}
	}
	
};

#endif