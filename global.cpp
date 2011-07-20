#include "global.h"

size_t nAutosomes = 22;
size_t nChromosomes = 24;

namespace name
{
	string common(const string& a, const string& b) {
		static const char specialChar = '.';
		
		size_t i, size = min(a.length(), b.length());
		for (i = 0; i < size; ++i) {
			if (a[i] != b[i]) break;
		}
		string s = "";
		if (i > 0) {
			// ignore special character
			if (i > 1 && a[i-1] == specialChar) --i;
			s = a.substr(0, i);
		}
		return s;
	}
	
	string fileext(const string& s) {
		size_t start = s.find_last_of('.');
		return s.substr(
			(start == string::npos) ? 0 : start + 1
		);
	}
	
	string filestem(const string& s) {
		size_t start = s.find_last_of('/');
		start = (start == string::npos) ? 0 : (start + 1);
		size_t end = s.find_last_of('.');
		return s.substr(
			start, 
			(end == string::npos) ? (s.length() - start + 1) : (end - start)
		);
	}
}