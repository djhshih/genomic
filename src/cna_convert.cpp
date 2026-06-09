#include "cna_convert.hpp"

Convert::Convert()
: Command("convert file format") {
	// Delcare options
	opts.add_options()
		("help", "print help message")
		("input,i", po::value< std::vector<std::string> >(), "input files")
		("from,f", po::value<std::string>(), "input file format [default: determined from file extension]")
		("output,o", po::value<std::string>(), "output file")
		("to,t", po::value<std::string>(), "output file format [default: either seg(as) or cn(as), depending on input file format]")
		;
	popts.add("input", -1);
}

void Convert::run() {
	if (vm.count("help")) {
		std::cout << "usage:  " << progname << " convert [options] <input files>" << std::endl;
		std::cout << opts << std::endl;
		return;
	}
	
	getOptions();
	
	switch (inputType) {
		case cna::data::segmented: {
			cna::SegmentedSampleSet<rvalue> set;
			if (outputType != cna::data::segmented_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case cna::data::segmented:
					set.write(outputFileName);
					break;
				case cna::data::segmented_ref: {
					cna::ReferenceSegmentedSampleSet<rvalue> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case cna::data::raw: {
					cna::RawSampleSet<rvalue> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw cna::invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case cna::data::segmented_ascn: {
			cna::SegmentedSampleSet<alleles_cn> set;
			if (outputType != cna::data::segmented_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case cna::data::segmented_ascn:
					set.write(outputFileName);
					break;
				case cna::data::segmented_ref: {
					cna::ReferenceSegmentedSampleSet<alleles_cn> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case cna::data::raw: {
					cna::RawSampleSet<alleles_cn> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw cna::invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case cna::data::raw: {
			cna::RawSampleSet<rvalue> set;
			if (outputType != cna::data::raw_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case cna::data::raw:
					set.write(outputFileName);
					break;
				case cna::data::raw_ref: {
					cna::ReferenceRawSampleSet<rvalue> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case cna::data::segmented: {
					cna::SegmentedSampleSet<rvalue> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw cna::invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case cna::data::raw_ascn: {
			cna::RawSampleSet<alleles_cn> set;
			if (outputType != cna::data::raw_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case cna::data::raw:
					set.write(outputFileName);
					break;
				case cna::data::raw_ref: {
					cna::ReferenceRawSampleSet<alleles_cn> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case cna::data::segmented: {
					cna::SegmentedSampleSet<alleles_cn> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw cna::invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		default:
			throw std::logic_error("Unsupported input type in convert.");
	}
}

void Convert::getOptions() {
	if (vm.count("input")) {
		inputFileNames = vm["input"].as< std::vector<std::string> >();
	} else {
		throw std::invalid_argument("No input file specified.");
	}
	
	if (vm.count("from")) {
		inputType = cna::mapping::extension[ vm["from"].as<std::string>() ];
	} else {
		inputType = cna::mapping::extension[ cna::name::fileext(inputFileNames[0]) ];
	}
	if (inputType == cna::data::invalid) {
		throw std::invalid_argument("Invalid input format type: '" + cna::name::fileext(inputFileNames[0]) + "'.");
	}
	
	cna::data::Type defaultOutputType = cna::data::invalid;
	switch (inputType) {
		case cna::data::raw:
		case cna::data::raw_ascn:
			defaultOutputType = cna::data::segmented_ascn;
			break;
		case cna::data::segmented:
			defaultOutputType = cna::data::raw;
			break;
		case cna::data::segmented_ascn:
			defaultOutputType = cna::data::raw_ascn;
			break;
		default:
			break;
	}
	
	if (vm.count("output")) {
		outputFileName = vm["output"].as<std::string>();
	} else {
		outputFileName = cna::name::filestem(inputFileNames[0]);
		if (vm.count("to")) {
			// use specified output format extension
			outputFileName += "." + vm["to"].as<std::string>();
		} else {
			// use default output format
			if (defaultOutputType != cna::data::invalid) {
				outputFileName += "." + cna::mapping::extension[defaultOutputType];
			}
		}
	}
	
	if (vm.count("to")) {
		outputType = cna::mapping::extension[ vm["to"].as<std::string>() ];
	} else {
		// determine output type from extension
		// outputFileName must be previously defined
		outputType = cna::mapping::extension[ cna::name::fileext(outputFileName) ];
	}
	if (outputType == cna::data::invalid) {
		throw std::invalid_argument("Invalid output format type for output file '" + outputFileName + "'.");
	}
}
