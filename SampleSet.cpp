#include "SampleSet.hpp"

void SampleSet::read(const vector<string>& fileNames, const string& markersFileName, const string& platform, bool isSorted) {
	if (fileNames.size() < 1) return;
	
	// Read markers, do not sort markers yet
	markers = marker::manager.create(platform);
	if (!markersFileName.empty()) {
		markers->setIO(delim, headerLine, nSkippedLines);
		markers->read(markersFileName, platform, false);
	}
	
	vector<string>::const_iterator it, end = fileNames.end();
	for (it = fileNames.begin(); it < end; ++it) {
		// read samples, append to set
		read(*it, platform, true);
	}
	
	// Sort after all samples have been read
	if (!markersFileName.empty()) markers->distribute();
	if (!isSorted) sort();
}

void SampleSet::read(const string& fileName, const string& platform, bool append) {
	file.open(fileName.c_str(), ios::in);
	if (!file.is_open()) throw runtime_error("Failed to open input file.");
	read(file, platform, fileName, append);
	file.close();
	trace("Read file %s\n", fileName.c_str());
}

void SampleSet::read(fstream& file, const string& platform, const string& fileName, bool append) {
	if (!append) clear();
	markers = marker::manager.create(platform);
	this->fileName = fileName;
	
	_read(file);
	if (!append) sort();
}

void SampleSet::write(const string& fileName) {
	file.open(fileName.c_str(), ios::out);
	if (!file.is_open()) throw runtime_error("Failed to open output file");
	write(file, fileName);
	file.close();
	trace("Wrote file %s\n", fileName.c_str());
}
