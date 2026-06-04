#ifndef genomic_sort_h
#define genomic_sort_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "genomic_common.hpp"
#include "Marker.hpp"


class Sort : public Command {
	
public:
	
	Sort() 
	: Command("sort markers file by genomic position") {
		
		// Delcare options
		opts.add_options()
			("help", "print help message")
			("input,i", po::value<string>(), "markers file")
			("output,o", po::value<string>(), "output sorted file")
			("named,m", po::value<bool>(), "whether marker names are present")
			;
		popts.add("input", 1).add("output", 1);
		
	}
	
	void run() {
		
		if (vm.count("help")) {
			cout << "usage:  " << progname << " sort [options] <markers file> <output file>" << endl;
			cout << opts << endl;
			return;
		}
		
		getOptions();
		
		string platform = "";
		
		marker::Set set(platform);
		
		// read and sort markers
		// all recognized chromosomes are discarded
		set.read(inputFileName, platform, true, named);
		
		// write markers to file
		set.write(outputFileName);
		
	}
	
private:
	
	string inputFileName, outputFileName;
	bool named;
	
	void getOptions() {
		
		if (vm.count("input")) {
			inputFileName = vm["input"].as<string>();
		} else {
			throw invalid_argument("Input file not specified.");
		}
		if (vm.count("output")) {
			outputFileName = vm["output"].as<string>();
		} else {
			outputFileName = name::filestem(inputFileName) + ".sorted." + name::fileext(inputFileName);
		}
		
		if (vm.count("named")) {
			named = vm["named"].as<bool>();
		} else {
			named = true;
		}
		
	}
	
};

#endif
