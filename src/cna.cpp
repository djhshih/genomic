#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "cna.hpp"

int main(int argc, char **argv)
{
	typedef map<string, reference_wrapper<Command> > CommandMap;
	
	Command* command = NULL;
	Convert convert;
	Filter filter;
	Clean clean;
	Sort sort;
	CommandMap commands;
	
	bool printUsage = false;
	bool printVersion = false;
	int ret = 0;
	
	try {
	
		commands.emplace("convert", ref(convert));
		commands.emplace("filter", ref(filter));
		commands.emplace("clean", ref(clean));
		commands.emplace("sort", ref(sort));
		
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
			
				CommandMap::const_iterator it = commands.find(s);
				if (it != commands.end()) {
					command = &(it->second.get());
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
			
			CommandMap::const_iterator it;
			CommandMap::const_iterator end = commands.end();
			for (it = commands.begin(); it != end; ++it) {
				cout << "  " << it->first << "  \t"
					<< it->second.get() << endl;
			}
			cout << endl;
			
			cout << "See '" << progname << " <command> --help' for more information on a specific command." << endl;
			
		}
		
		if (printVersion) {
			cout << "cna version "
				<< cna_VERSION_MAJOR << '.' << cna_VERSION_MINOR
				<< endl;
		}
	
	} catch (const exception& e) {
		cout << e.what() << endl;
		ret = 1;
	}
	
	return ret;
}
