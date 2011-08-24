#ifndef genomic_FilesDiff_h
#define genomic_FilesDiff_h

#include <fstream>
#include <string>
#include <stdexcept>

using namespace std;

class FilesDiff
{
private:
	fstream file1, file2;
	size_t nSkippedLines;
	
public:
	FilesDiff() : nSkippedLines(0) {}
	FilesDiff(size_t skip) : nSkippedLines(skip) {}
	~FilesDiff() {
		if (file1.is_open()) file2.close();
		if (file1.is_open()) file2.close();
	}
	int different(const string&fileName1, const string& fileName2) {
		file1.open(fileName1.c_str(), ios::in);
		file2.open(fileName2.c_str(), ios::in);
		if (!file1.is_open() || !file2.is_open()) {
			throw runtime_error("Faile to open files for comparison");
		}
		
		string s1, s2;
		int diff = 0;
		size_t lineCount = 0;
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