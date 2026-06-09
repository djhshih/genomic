#include "cna_filter.hpp"

Filter::Filter()
: Command("filter file based on reference") {
	// Delcare options
	opts.add_options()
		("help", "print help message")
		("input,i", po::value<std::string>(), "sample set file")
		("reference,r", po::value<std::string>(), "reference set file (same format as sample set)")
		("format,f", po::value<std::string>(), "input file format [default: determined from file extension]")
		("reference_format,g", po::value<std::string>(), "reference file format [default: reference format corresponding to input]")
		("output,o", po::value<std::string>(), "output file")
		("threshold,t", po::value<float>(), "overlap threshold (only used for segmentation files) [default: 0.5]")
		("score,s", po::value<std::string>(), "overlap score method, only used for segmentation files) [options: dice (default), query, reference, min, max]")
		("inverse,v", po::value<bool>(), "select instead of filter overlapping segments")
		("merge,m", po::value<bool>(), "merge filtered segments with upstream/downstream segments?")
		("aberrant,a", po::value<bool>(), "filter only aberrant segments?")
		("state_diff", po::value<rvalue>(), "threshold for difference from reference state")
		("ref_state", po::value<rvalue>(), "reference state")
		("optimize,O", po::value<bool>(), "optimize algorithm speed, assuming contiguity of reference segments")
		;
	popts.add("input", 1).add("reference", 1).add("output", 1);
}

void Filter::run() {
	if (vm.count("help")) {
		std::cout << "usage:  " << progname << " filter [options] <sample set file> <reference set file> <output file>" << std::endl;
		std::cout << opts << std::endl;
		return;
	}
	
	getOptions();
	
	switch (inputType) {
		case cna::data::segmented:
			switch(referenceType) {
				case cna::data::segmented_ref:
					filter< cna::SegmentedSampleSet<rvalue>, cna::ReferenceSegmentedSampleSet<rvalue> >(segmented<true>());
					break;
				case cna::data::segmented:
					filter< cna::SegmentedSampleSet<rvalue>, cna::SegmentedSampleSet<rvalue> >(segmented<true>());
					break;
				default:
					throw std::invalid_argument("Unsupported reference file format.");
			}
			break;
		case cna::data::segmented_ascn:
			switch(referenceType) {
				case cna::data::segmented_ref:
					filter< cna::SegmentedSampleSet<alleles_cn>, cna::ReferenceSegmentedSampleSet<alleles_cn> >(segmented<true>());
					break;
				case cna::data::segmented_ascn:
					filter< cna::SegmentedSampleSet<alleles_cn>, cna::SegmentedSampleSet<alleles_cn> >(segmented<true>());
					break;
				default:
					throw std::invalid_argument("Unsupported reference file format.");
			}
			break;
		case cna::data::raw:
			switch(referenceType) {
				case cna::data::raw_ref:
					filter< cna::RawSampleSet<rvalue>, cna::ReferenceRawSampleSet<rvalue> >(segmented<false>());
					break;
				case cna::data::raw:
					filter< cna::RawSampleSet<rvalue>, cna::RawSampleSet<rvalue> >(segmented<false>());
					break;
				default:
					throw std::invalid_argument("Unsupported reference file format.");
			}
			break;
		case cna::data::raw_ascn:
			switch(referenceType) {
				case cna::data::raw_ref:
					filter< cna::RawSampleSet<alleles_cn>, cna::ReferenceRawSampleSet<alleles_cn> >(segmented<false>());
					break;
				case cna::data::raw_ascn:
					filter< cna::RawSampleSet<alleles_cn>, cna::RawSampleSet<alleles_cn> >(segmented<false>());
					break;
				default:
					throw std::invalid_argument("Unsupported reference file format.");
			}
			break;
		default:
			throw std::logic_error("Unsupported input type in filter.");
	}
}

void Filter::getOptions() {
	if (vm.count("input") && vm.count("reference")) {
		inputFileName = vm["input"].as<std::string>();
		referenceFileName = vm["reference"].as<std::string>();
	} else {
		throw std::invalid_argument("Both input and reference files must be specified.");
	}
	
	if (vm.count("format")) {
		inputType = cna::mapping::extension[ vm["format"].as<std::string>() ];
	} else {
		inputType = cna::mapping::extension[ cna::name::fileext(inputFileName) ];
	}
	if (inputType == cna::data::invalid) {
		throw std::invalid_argument("Invalid input format type for input file '" + inputFileName + "'.");
	}
	
	if (vm.count("reference_format")) {
		referenceType = cna::mapping::extension[ vm["reference_format"].as<std::string>() ];
	} else {
		if (inputType == cna::data::raw || inputType == cna::data::raw_ascn || inputType == cna::data::raw_lrrbaf) {
			referenceType = cna::data::raw_ref;
		} else {
			referenceType = cna::data::segmented_ref;
		}
	}
	if (referenceType == cna::data::invalid) {
		throw std::invalid_argument("Invalid reference format type for reference file '" + referenceFileName + "'.");
	}
	
	if (vm.count("output")) {
		outputFileName = vm["output"].as<std::string>();
	} else {
		outputFileName = cna::name::filestem(inputFileName) + ".filtered." + cna::name::fileext(inputFileName);
	}
	
	if (vm.count("threshold")) {
		threshold = vm["threshold"].as<float>();
	} else {
		threshold = 0.5;
	}
	
	if (vm.count("score")) {
		score = vm["score"].as<std::string>();
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
