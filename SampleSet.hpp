#ifndef genomic_SampleSet_h
#define genomic_SampleSet_h

#include <string>
#include <stdexcept>
#include <fstream>

#include "global.hpp"
#include "Sample.hpp"
#include "Marker.hpp"

class GenericSampleSet;

class SampleSet
{
	friend class GenericSampleSet;
	//  for accessing the private clone() function
	//    and the read() and write() functions
public:
	SampleSet() : markers(NULL) {
		setIO();
	}
	SampleSet(marker::Set* markerSet) : markers(markerSet) {
		setIO();
	}
	virtual ~SampleSet() {
		if (file.is_open()) file.close();
	}
	void setIO(char _delim='\t', size_t _headerLine=1, size_t _nSkippedLines=0) {
		delim = _delim;
		headerLine = _headerLine;
		nSkippedLines = _nSkippedLines;
	};
	void read(const vector<string>& fileNames, const string& markersFileName, bool isSorted=false) {
		read(fileNames, markersFileName, markersFileName, isSorted);
	}
	void read(const vector<string>& fileNames, const string& markersFileName, const string& platform, bool isSorted) {
		if (fileNames.size() < 1) return;
		
		// Read markers, do not sort markers yet
		markers = marker::manager.create(platform);
		markers->setIO(delim, headerLine, nSkippedLines);
		markers->read(markersFileName, platform, false);
		
		vector<string>::const_iterator it, end = fileNames.end();
		for (it = fileNames.begin(); it < end; ++it) {
			// read samples, append to set
			read(*it, platform, true);
		}
		
		// Sort after all samples have been read
		markers->distribute();
		if (!isSorted) sort();
	}
	void read(const string& fileName, bool append=false) {
		// use fileName also as platform name
		read(fileName, fileName, append);
	}
	void read(const string& fileName, const string& platform, bool append=false) {
		file.open(fileName.c_str(), ios::in);
		if (!file.is_open()) throw runtime_error("Failed to open input file");
		read(file, platform, fileName, append);
		file.close();
		trace("Read file %s\n", fileName.c_str());
	}
	void read(fstream& file, const string& platform, const string& fileName, bool append=false) {
		if (!append) clear();
		markers = marker::manager.create(platform);
		//this->fileName = "." + fileType;
		this->fileName = fileName;
		
		_read(file);
		if (!append) sort();
	}
	void write(const string& fileName) {
		file.open(fileName.c_str(), ios::out);
		if (!file.is_open()) throw runtime_error("Failed to open output file");
		write(file, fileName);
		file.close();
		trace("Wrote file %s\n", fileName.c_str());
	}
	void write(fstream& file, const string& fileName) {
		this->fileName = fileName;
		_write(file);
	}
	virtual void clear() = 0;
	virtual void sort() = 0;
	virtual data::Type type() = 0;
	virtual size_t size() = 0;
	
protected:
	char delim;
	size_t nSkippedLines;
	size_t headerLine;
	bool mergeSamples;
	
	string fileName;
	marker::Set* markers;
private:
	fstream file;
	//  file IO is the responsibiliity of the base class
	virtual void _read(fstream& file) = 0;
	virtual void _write(fstream& file) = 0;
	virtual SampleSet* clone() = 0;
};

#endif