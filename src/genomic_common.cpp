#include "genomic_common.hpp"

std::string progname = "genomic";

std::ostream& operator<<(std::ostream& os, const Command& c) {
	return os << c.description;
}
