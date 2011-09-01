#ifndef genomic_fiilter_h
#define genomic_fiilter_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "genomic_common.hpp"


class Convert : public Command {
public:
	Convert()
	: Command("convert file format") {
		// Delcare options
		//po::options_description opts("Convert options");
		opts.add_options()
			("help", "print help message")
			("input,i", po::value< vector<string> >(), "input files")
			("from,f", po::value<string>(), "input file format [default: determined from file extension]")
			("output,o", po::value<string>(), "output file")
			("to,t", po::value<string>(), "output file format [default: either seg(as) or cn(as), depending on input file format]")
			;
		popts.add("input", -1);
	}
	
	void run() {
		
		if (vm.count("help")) {
			cout << "usage:  " << progname << " convert [options] <input files>" << endl;
			cout << opts << endl;
			return;
		}
		
		getOptions();
		
		
		switch (inputType) {
			
			case data::picnic: {
				
				const data::Type defaultOutType = data::segmented_ascn;
		
				const string& markersFileName = inputFileNames[0];
				const string& samplesFileName = inputFileNames[1];
				
				cout << "Input files: " << inputFileNames << endl;

				ifstream samplesFile(samplesFileName.c_str());
				if (!samplesFile.is_open()) {
					throw runtime_error("Failed to open samples file.");
				}
				
				vector<string> fileNames;
				while (!samplesFile.eof()) {
					string sample;
					samplesFile >> sample;
					if (sample != "") {
						fileNames.push_back(name::filepath(samplesFileName) + sample);
						cout << sample << endl;
					}
				}
					
				PicnicSampleSet pset;
				pset.read(fileNames, markersFileName, true);
				
				switch (outputType) {
					case data::segmented:
					case data::segmented_ascn: {
						SegmentedSampleSet<PicnicSampleSet::Value> sset(pset);
						sset.write(outputFileName);
						break;
					}
					case data::raw:
					case data::raw_ascn: {
						pset.write(outputFileName);
						break;
					}
					default:
						throw invalid_conversion("Unsupported output type.");
				}
				
				break;
			}
			
			case data::dchip: {
				
				throw logic_error("Conversion not implemented");
				break;
			}
			
			case data::penncnv: {
				
				throw logic_error("Conversion not implemented");
				break;
			}
				
			case data::cnag: {
					
				throw logic_error("Conversion not implemented");
				break;
			}
					
			case data::segmented: {
			
				SegmentedSampleSet<rvalue> set;
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::segmented:
						set.write(outputFileName);
						break;
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
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::segmented_ascn:
						set.write(outputFileName);
						break;
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
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::raw:
						set.write(outputFileName);
						break;
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
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::raw:
						set.write(outputFileName);
						break;
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
					
		}
		
	}
	
private:
	
	vector<string> inputFileNames;
	string outputFileName;
	data::Type inputType, outputType;

	void getOptions() {
		
		if (vm.count("input")) {
			inputFileNames = vm["input"].as< vector<string> >();
		} else {
			throw invalid_argument("No input file specified.");
		}
		
		if (vm.count("from")) {
			inputType = mapping::extension[ vm["from"].as<string>() ];
		} else {
			inputType = mapping::extension[ name::fileext(inputFileNames[0]) ];
		}
		if (inputType == data::invalid) {
			throw invalid_argument("Invalid input format type.");
		}
		
		data::Type defaultOutputType = data::invalid;
		switch (inputType) {
			case data::raw:
			case data::cnag:
			case data::dchip:
			case data::penncnv:
				defaultOutputType = data::segmented;
				break;
			case data::picnic:
			case data::raw_ascn:
				defaultOutputType = data::segmented_ascn;
				break;
			case data::segmented:
				defaultOutputType = data::raw;
				break;
			case data::segmented_ascn:
				defaultOutputType = data::raw_ascn;
				break;
		}
		
		if (vm.count("output")) {
			outputFileName = vm["output"].as<string>();
		} else {
			outputFileName = name::filestem(inputFileNames[0]);
			if (defaultOutputType != data::invalid) {
				outputFileName += "." + mapping::extension[defaultOutputType];
			}
		}
		
		if (vm.count("to")) {
			outputType = mapping::extension[ vm["to"].as<string>() ];
		} else {
			outputType = mapping::extension[ name::fileext(outputFileName) ];
		}
		if (outputType == data::invalid) {
			throw invalid_argument("Invalid output format type.");
		}
	}
	
};

#endif
