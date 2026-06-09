#include "cna_common.hpp"

std::string progname = "cna";

std::ostream& operator<<(std::ostream& os, const Command& c) {
	return os << c.description;
}
