#ifndef genomic_convert_h
#define genomic_convert_h

#include <stdexcept>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "genomic_common.hpp"


class Convert : public Command {
public:
	Convert();
	
	void run();
	
private:
	
	std::vector<std::string> inputFileNames;
	std::string outputFileName;
	data::Type inputType, outputType;

	void getOptions();
	
};

#endif
