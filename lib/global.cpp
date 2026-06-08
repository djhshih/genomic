#include "global.hpp"

namespace cna {

chromid nAutosomes = 22;
chromid nChromosomes = 24;

namespace mapping {
	ChromosomesMap chromosome;
	ExtensionMap extension;
}

namespace name
{
	std::string common(const std::string& a, const std::string& b) {
		static const char specialChar = '.';
		
		size_t i, size = std::min(a.length(), b.length());
		for (i = 0; i < size; ++i) {
			if (a[i] != b[i]) break;
		}
		std::string s = "";
		if (i > 0) {
			// ignore special character
			if (i > 1 && a[i-1] == specialChar) --i;
			s = a.substr(0, i);
		}
		return s;
	}
	
	std::string fileext(const std::string& s) {
		size_t start = s.find_last_of('.');
		return s.substr(
			(start == std::string::npos || start == s.size()-1) ? 0 : start + 1
		);
	}
	
	std::string filestem(const std::string& s) {
		size_t start = s.find_last_of('/');
		start = (start == std::string::npos) ? 0 : (start + 1);
		size_t end = s.find_last_of('.');
		return s.substr(
			start, 
			(end == std::string::npos) ? (s.length() - start + 1) : (end - start)
		);
	}
	
	std::string filepath(const std::string& s) {
		size_t start = 0;
		size_t end = s.find_last_of('/');
		if (end == std::string::npos) {
			return "./";
		} else {
			return s.substr(start, end+1);
		}
	}
		
}

} // namespace cna
