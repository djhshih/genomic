#include "parse.hpp"

#include <cctype>

FieldScanner::FieldScanner(const std::string& text, char delimiter)
	: line(text), delim(delimiter), pos(0), whitespaceMode(delimiter == ' ') {}

bool FieldScanner::next(std::string_view& field) {
	if (whitespaceMode) {
		while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos])) != 0) {
			++pos;
		}
		if (pos >= line.size()) {
			return false;
		}
		const size_t start = pos;
		while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos])) == 0) {
			++pos;
		}
		field = line.substr(start, pos - start);
		return true;
	}

	if (pos > line.size()) {
		return false;
	}
	const size_t start = pos;
	while (pos < line.size() && line[pos] != delim) {
		++pos;
	}
	field = line.substr(start, pos - start);
	if (pos < line.size()) {
		++pos;
	} else {
		pos = line.size() + 1;
	}
	return true;
}
