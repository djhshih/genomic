#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "genomic.hpp"

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
