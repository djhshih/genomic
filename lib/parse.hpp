#ifndef genomic_parse_h
#define genomic_parse_h

#include <string>
#include <string_view>
#include <charconv>

class FieldScanner {
private:
	std::string_view line;
	char delim;
	size_t pos;
	bool whitespaceMode;

public:
	FieldScanner(const std::string& text, char delimiter);
	bool next(std::string_view& field);
};

template <typename T> inline
bool parseNumber(std::string_view text, T& value) {
	const char* begin = text.data();
	const char* end = begin + text.size();
	const std::from_chars_result result = std::from_chars(begin, end, value);
	return result.ec == std::errc() && result.ptr == end;
}

#endif
