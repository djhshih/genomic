#ifndef genomic_command_h
#define genomic_command_h

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

class missing_optional : public std::runtime_error {
public:
	explicit missing_optional(const std::string& what_arg)
	: std::runtime_error("Unspecified optional: " + what_arg) {}
};

class missing_argument : public std::runtime_error {
public:
	explicit missing_argument(const std::string& what_arg)
	: std::runtime_error("Missing argument: " + what_arg) {}
};

extern std::string progname;

class Command {
public:
	Command() {}
	Command(const std::string& desc) : description(desc) {}
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
	std::string description;
	
protected:
	po::variables_map vm;
	po::options_description opts;
	po::positional_options_description popts;
};


std::ostream& operator<<(std::ostream& os, const Command& c);

// helper function to print vectors
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}

#endif
