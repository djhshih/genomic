#include "cna_clean.hpp"

Clean::Clean()
: Command("clean file based on filtering criteria") {
	// Delcare options
	opts.add_options()
		("help", "print help message")
		("input,i", po::value<std::string>(), "sample set file")
		("output,o", po::value<std::string>(), "output file")
		("format,f", po::value<std::string>(), "input file format [default: determined from file extension]")
		("inverse,v", po::value<bool>(), "select instead of filter overlapping segments")
		("merge,m", po::value<bool>(), "merge filtered segments with upstream/downstream segments?")
		("count", po::value<position>(), "threshold for number of markers in segment")
		("length", po::value<position>(), "threshold for segment length")
		("balanced", po::value<bool>(), "remove balanced segments?")
		("state_diff", po::value<rvalue>(), "threshold for difference from reference state")
		("ref_state", po::value<rvalue>(), "reference state")
		;
	popts.add("input", 1).add("output", 1);
}

void Clean::run() {
	if (vm.count("help")) {
		std::cout << "usage:  " << progname << " clean [options] <sample set file> <output file>" << std::endl;
		std::cout << opts << std::endl;
		return;
	}
	
	getOptions();
	
	switch (inputType) {
		case cna::data::segmented:
			clean< cna::SegmentedSampleSet<rvalue> >(segmented<true>());
			break;
		case cna::data::segmented_ascn:
			clean< cna::SegmentedSampleSet<alleles_cn> >(segmented<true>());
			break;
		case cna::data::raw:
			clean< cna::RawSampleSet<rvalue> >(segmented<false>());
			break;
		case cna::data::raw_ascn:
			clean< cna::RawSampleSet<alleles_cn> >(segmented<false>());
			break;
		default:
			throw std::logic_error("Unsupported input type in clean.");
	}
}

void Clean::getOptions() {
	if (vm.count("input")) {
		inputFileName = vm["input"].as<std::string>();
	} else {
		throw std::invalid_argument("Input file not specified for clean command.");
	}
	
	if (vm.count("format")) {
		inputType = cna::mapping::extension[ vm["format"].as<std::string>() ];
	} else {
		inputType = cna::mapping::extension[ cna::name::fileext(inputFileName) ];
	}
	if (inputType == cna::data::invalid) {
		throw std::invalid_argument("Invalid input format type for input file '" + inputFileName + "'.");
	}
	
	if (vm.count("output")) {
		outputFileName = vm["output"].as<std::string>();
	} else {
		outputFileName = cna::name::filestem(inputFileName) + ".filtered." + cna::name::fileext(inputFileName);
	}

	if (vm.count("inverse")) {
		inverse = vm["inverse"].as<bool>();
	} else {
		inverse = false;
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
