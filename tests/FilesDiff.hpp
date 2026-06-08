#ifndef cna_FilesDiff_h
#define cna_FilesDiff_h

#include <fstream>
#include <string>
#include <stdexcept>

class FilesDiff
{
private:
	std::fstream file1, file2;
	std::size_t nSkippedLines;
	
public:
	FilesDiff() : nSkippedLines(0) {}
	FilesDiff(std::size_t skip) : nSkippedLines(skip) {}
	~FilesDiff() {
		if (file1.is_open()) file2.close();
		if (file1.is_open()) file2.close();
	}
	int different(const std::string&fileName1, const std::string& fileName2) {
		file1.open(fileName1.c_str(), std::ios::in);
		file2.open(fileName2.c_str(), std::ios::in);
		if (!file1.is_open() || !file2.is_open()) {
			throw std::runtime_error("Failed to open files for comparison: '" + fileName1 + "' and '" + fileName2 + "'.");
		}
		
		std::string s1, s2;
		int diff = 0;
		std::size_t lineCount = 0;
		while(true) {
			getline(file1, s1);
			getline(file2, s2);
			if (file1.eof()) {
				if (!file2.eof()) {
					// file1 is shorter
					diff = -1;  
				}
				break;
			} else if (file2.eof()) {
				// file2 is shorter
				diff = -2;
				break;
			}
			// do not check for different during line skipping
			if (++lineCount <= nSkippedLines) continue;
			if (s1 != s2) {
				diff = 1;
				break;
			}
		}
		
		file1.close();
		file2.close();
		
		return diff;
	}
};	

#endif