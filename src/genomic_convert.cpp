#include "genomic_convert.hpp"

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
		case data::segmented: {
			SegmentedSampleSet<rvalue> set;
			if (outputType != data::segmented_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case data::segmented:
					set.write(outputFileName);
					break;
				case data::segmented_ref: {
					ReferenceSegmentedSampleSet<rvalue> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case data::raw: {
					RawSampleSet<rvalue> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case data::segmented_ascn: {
			SegmentedSampleSet<alleles_cn> set;
			if (outputType != data::segmented_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case data::segmented_ascn:
					set.write(outputFileName);
					break;
				case data::segmented_ref: {
					ReferenceSegmentedSampleSet<alleles_cn> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case data::raw: {
					RawSampleSet<alleles_cn> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case data::raw: {
			RawSampleSet<rvalue> set;
			if (outputType != data::raw_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case data::raw:
					set.write(outputFileName);
					break;
				case data::raw_ref: {
					ReferenceRawSampleSet<rvalue> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case data::segmented: {
					SegmentedSampleSet<rvalue> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw invalid_conversion("Unsupported output type.");
			}
			
			break;
		}
		case data::raw_ascn: {
			RawSampleSet<alleles_cn> set;
			if (outputType != data::raw_ref) {
				set.read(inputFileNames);
			}
			
			switch (outputType) {
				case data::raw:
					set.write(outputFileName);
					break;
				case data::raw_ref: {
					ReferenceRawSampleSet<alleles_cn> out;
					out.read(inputFileNames);
					out.write(outputFileName);
					break;
				}
				case data::segmented: {
					SegmentedSampleSet<alleles_cn> out(set);
					out.write(outputFileName);
					break;
				}
				default:
					throw invalid_conversion("Unsupported output type.");
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
		inputType = mapping::extension[ vm["from"].as<std::string>() ];
	} else {
		inputType = mapping::extension[ name::fileext(inputFileNames[0]) ];
	}
	if (inputType == data::invalid) {
		throw std::invalid_argument("Invalid input format type: '" + name::fileext(inputFileNames[0]) + "'.");
	}
	
	data::Type defaultOutputType = data::invalid;
	switch (inputType) {
		case data::raw:
		case data::raw_ascn:
			defaultOutputType = data::segmented_ascn;
			break;
		case data::segmented:
			defaultOutputType = data::raw;
			break;
		case data::segmented_ascn:
			defaultOutputType = data::raw_ascn;
			break;
		default:
			break;
	}
	
	if (vm.count("output")) {
		outputFileName = vm["output"].as<std::string>();
	} else {
		outputFileName = name::filestem(inputFileNames[0]);
		if (vm.count("to")) {
			// use specified output format extension
			outputFileName += "." + vm["to"].as<std::string>();
		} else {
			// use default output format
			if (defaultOutputType != data::invalid) {
				outputFileName += "." + mapping::extension[defaultOutputType];
			}
		}
	}
	
	if (vm.count("to")) {
		outputType = mapping::extension[ vm["to"].as<std::string>() ];
	} else {
		// determine output type from extension
		// outputFileName must be previously defined
		outputType = mapping::extension[ name::fileext(outputFileName) ];
	}
	if (outputType == data::invalid) {
		throw std::invalid_argument("Invalid output format type for output file '" + outputFileName + "'.");
	}
}
