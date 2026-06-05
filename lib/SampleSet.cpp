#include "SampleSet.hpp"

void SampleSet::read(const std::vector<std::string>& fileNames, const std::string& markersFileName, const std::string& platform, bool isSorted) {
	if (fileNames.size() < 1) return;
	
	// Read markers, do not sort markers yet
	markers = marker::manager.create(platform);
	if (!markersFileName.empty()) {
		markers->setIO(io);
		markers->read(markersFileName, platform, false);
	}
	
	std::vector<std::string>::const_iterator it, end = fileNames.end();
	for (it = fileNames.begin(); it < end; ++it) {
		// read samples, append to set
		read(*it, platform, true);
	}
	
	// Sort after all samples have been read
	if (!markersFileName.empty()) markers->distribute();
	if (!isSorted) sort();
}

void SampleSet::read(const std::string& fileName, const std::string& platform, bool append) {
	file.open(fileName.c_str(), std::ios::in);
	if (!file.is_open()) throw std::runtime_error("Failed to open input file '" + fileName + "'.");
	read(file, platform, fileName, append);
	file.close();
	log_trace(__FILE__, __LINE__, __func__, "Read file %s", fileName.c_str());
}

void SampleSet::read(std::fstream& file, const std::string& platform, const std::string& fileName, bool append) {
	if (!append) clear();
	markers = marker::manager.create(platform);
	this->fileName = fileName;
	
	_read(file);
	if (!append) sort();
}

void SampleSet::write(const std::string& fileName) {
	file.open(fileName.c_str(), std::ios::out);
	if (!file.is_open()) throw std::runtime_error("Failed to open output file '" + fileName + "'.");
	write(file, fileName);
	file.close();
	log_trace(__FILE__, __LINE__, __func__, "Wrote file %s", fileName.c_str());
}
