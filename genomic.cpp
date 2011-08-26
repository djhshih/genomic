#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "genomic.hpp"


string progname = "genomic";

// helper function to print vectors
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}

class missing_optional : public std::runtime_error {
public:
	explicit missing_optional(const string& what_arg)
	: std::runtime_error("Unspecified optional: " + what_arg) {}
};

class missing_argument : public std::runtime_error {
public:
	explicit missing_argument(const string& what_arg)
	: std::runtime_error("Missing argument: " + what_arg) {}
};

class Command {
public:
	Command() {}
	Command(const string& desc) : description(desc) {}
	virtual ~Command() {}
	
	Command& parse(int argc, char *argv[]) {
		store( po::command_line_parser(argc, argv)
					.options(opts)
					.positional(popts)
					.run(), vm );
		notify(vm);
		return *this;
	}
	
	virtual void run() = 0;
	
public:
	string description;
	
protected:
	po::variables_map vm;
	po::options_description opts;
	po::positional_options_description popts;
};

ostream& operator<<(ostream& os, const Command& c) {
	return os << c.description;
}

class Convert : public Command {
public:
	Convert()
	: Command("convert file format") {
		// Delcare options
		//po::options_description opts("Convert options");
		opts.add_options()
			("help", "print help message")
			("input,i", po::value< vector<string> >(), "input files")
			("from,f", po::value<string>(), "input file format")
			("output,o", po::value<string>(), "output file")
			("to,t", po::value<string>(), "output file format")
			;
		popts.add("input", -1);
	}
	
	void run() {
		
		if (vm.count("help")) {
			cout << "usage:  " << progname << " convert <args>" << endl;
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
			
				SegmentedSampleSet<cnvalue> set;
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::segmented:
						set.write(outputFileName);
						break;
					case data::raw: {
						RawSampleSet<cnvalue> out(set);
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
				
				RawSampleSet<cnvalue> set;
				set.read(inputFileNames);
				
				switch (outputType) {
					case data::raw:
						set.write(outputFileName);
						break;
					case data::segmented: {
						SegmentedSampleSet<cnvalue> out(set);
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

	/*
	SampleSet* newSampleSet(data::Type type, SampleSet* set) {
		SampleSet* out;
		switch (type) {
			case data::segmented:
				out = new SegmentedSampleSet<cnvalue>(*set);
				break;
			case data::segmented_ascn:
				out = new SegmentedSampleSet<alleles_cn>(*set);
				break;
			case data::raw:
				out = new RawSampleSet<cnvalue>(*set);
				break;
			case data::raw_ascn:
				out = new RawSampleSet<alleles_cn>(*set);
				break;
		}
		return out;
	}
	*/
	
	void getOptions() {
		
		if (vm.count("input")) {
			inputFileNames = vm["input"].as< vector<string> >();
		} else {
			throw runtime_error("No input file specified.");
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

class Filter : public Command {
public:
	Filter() 
	: Command("filter file based on reference") {
		
		// Delcare options
		opts.add_options()
			("help", "print help message")
			("input,i", po::value< vector<string> >(), "input files")
			("output,o", po::value<string>(), "output file")
			;
		popts.add("input", 3).add("output", 1);
		
	}
	
	void run() {
		if (vm.count("help")) {
			cout << "usage:  " << progname << " filter <args>" << endl;
			cout << opts << endl;
			return;
		}
		
	}
	
};

int main(int argc, char **argv)
{
	typedef map<string, Command*> CommandMap;
	
	Command* command = NULL;
	CommandMap commands;
	
	bool printUsage = false;
	bool printVersion = false;
	int ret = 0;
	
	try {
	
		commands["convert"] = new Convert();
		commands["filter"] = new Filter();
		
		// Use the first argument (excluding name of program itself)
		//   to determine the command
		if (argc < 2) {
			printUsage = true;
		} else {
			string s = argv[1];
			
			if (s == "--help" || s == "-h") {
				printUsage = true;
			} else if (s == "--version" || s == "-v") {
				printVersion = true;
			} else {
			
				CommandMap::iterator it;
				CommandMap::const_iterator end = commands.end();
				for (it = commands.begin(); it != end; ++it) {
					if (s == it->first) {
						command = it->second;
						break;
					}
				}
				
				if (!command) {
					string err;
					err += "Error: '" + s + "' is an unrecognized command.";
					throw runtime_error(err);
				}
				
			}
			
		}
		
		if (command) {
			// forward remaining command line arguments
			//   for parsing by specified command module, and run the module
			(*command).parse(argc-1, argv+1).run();
		}
		
		if (printUsage) {
			cout << "usage:  " << progname << " <command> [<args>]" << endl;
			cout << endl;
			
			cout << "commands:" << endl;
			
			CommandMap::iterator it;
			CommandMap::const_iterator end = commands.end();
			for (it = commands.begin(); it != end; ++it) {
				cout << "  " << it->first << '\t'
					<< *(it->second) << endl;
			}
			cout << endl;
			
			cout << "See '" << progname << " <command> --help' for more information on a specific command." << endl;
			
		}
		
		if (printVersion) {
			cout << "genomic version "
				<< genomic_VERSION_MAJOR << '.' << genomic_VERSION_MINOR
				<< endl;
		}
	
	} catch (const exception& e) {
		cout << e.what() << endl;
		ret = 1;
	}
	
	// clear command objects
	CommandMap::iterator cmdIt;
	CommandMap::const_iterator cmdEnd = commands.end();
	for (cmdIt = commands.begin(); cmdIt != cmdEnd; ++cmdIt) {
		delete cmdIt->second;
	}

	return ret;
}
