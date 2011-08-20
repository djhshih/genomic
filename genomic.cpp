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

// helper function to print vectors
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}

string progname = "genomic";

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
		if (vm.count("input")) {
			const vector<string>& inputs = vm["input"].as< vector<string> >();
			cout << "Input files: " << inputs << endl;
			
			if (vm.count("from")) {
				
				if (vm["from"].as<string>() == "picnic") {
			
					const string& markersFileName = inputs[0];
					const string& samplesFileName = inputs[1];

					ifstream samplesFile(samplesFileName.c_str());
					if (!samplesFile.is_open()) {
						throw runtime_error("Failed to open samples file");
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
					//pset.write("output.picnic.ascn");
					
					string outFileName = "output.segas";
					if (vm.count("output")) {
						outFileName = vm["output"].as<string>();
					}
					
					SegmentedSampleSet<PicnicSampleSet::Value> sset(pset);
					sset.write(outFileName);
				}
				
				//TODO add other formats
				
			} else {
				throw runtime_error("Error: no input format specified");
			}
			
		} else {
			throw runtime_error("Error: no input file specified");
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
					err += "ERROR: '" + s + "' is unrecognized command";
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
	
	} catch (exception& e) {
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
