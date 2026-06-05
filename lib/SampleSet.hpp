#ifndef genomic_SampleSet_h
#define genomic_SampleSet_h

#include <string>
#include <stdexcept>
#include <fstream>

#include "global.hpp"
#include "Sample.hpp"
#include "Marker.hpp"
#include "Properties.hpp"

class GenericSampleSet;

//TODO use boost::spirit for parsing
//TODO allow use of input delimiter other than whitespace in all functions

class SampleSet
{
	friend class GenericSampleSet;
	//  for accessing the private clone() function
	//    and the read() and write() functions
	
public:
	
	SampleSet() : markers(NULL) {}
	
	SampleSet(marker::Set* markerSet) : markers(markerSet) {}
	
	SampleSet(const SampleSet& other)
	: io(other.io), markers(other.markers), fileName(other.fileName) {
		marker::manager.ref(markers);
	}
	
	virtual ~SampleSet() {
		if (file.is_open()) file.close();
	}
	
	void setIO(const IOProperties& io) {
		this->io = io;
	}
	
	void read(const std::vector<std::string>& fileNames, bool isSorted=false) {
		read(fileNames, "", isSorted);
	}
	
	void read(const std::vector<std::string>& fileNames, const std::string& markersFileName, bool isSorted=false) {
		std::string platform;
		marker::manager.newSetName(platform);
		read(fileNames, markersFileName, platform, isSorted);
	}
	
	void read(const std::vector<std::string>& fileNames, const std::string& markersFileName, const std::string& platform, bool isSorted);
	
	void read(const std::string& fileName, bool append=false) {
		std::string platform;
		marker::manager.newSetName(platform);
		read(fileName, platform, append);
	}
	
	// N.B. If platform name is specified explicitly by user, the marker file should be already sorted!
	void read(const std::string& fileName, const std::string& platform, bool append=false);
	
	void read(std::fstream& file, const std::string& platform, const std::string& fileName, bool append=false);
	
	void write(const std::string& fileName);
	
	void write(std::fstream& file, const std::string& fileName) {
		this->fileName = fileName;
		_write(file);
	}
	
	virtual void clear() = 0;
	virtual void sort() = 0;
	virtual data::Type type() = 0;
	virtual size_t size() = 0;
	
protected:
	
	IOProperties io;
	CNACriteria cna;
	
	std::string fileName;
	marker::Set* markers;
	
private:
	
	std::fstream file;
	//  file IO is the responsibiliity of the base class
	virtual void _read(std::fstream& file) = 0;
	virtual void _write(std::fstream& file) = 0;
	virtual SampleSet* clone() const = 0;
	
};

#endif
