#ifndef genomic_SplitRawSampleSet_h
#define genomic_SplitRawSampleSet_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>

#include "CopyNumberValue.hpp"
#include "SampleSet.hpp"

extern marker::Manager marker::manager;

template <typename V> class RawSampleSet;


template <typename V> 
class SplitRawSampleSet : public RawSampleSet<V>
{
public:
	typedef RawSampleSet<V> Base;
private:
	size_t dataColumn;
	void _read(fstream& file);
	void readSampleValue(istringstream& stream, typename Base::RawSample* sample, size_t chromIndex, const char delim);
public:
	SplitRawSampleSet() : dataColumn(1) {}
	SplitRawSampleSet(size_t sampleDataColumn) : dataColumn(sampleDataColumn) {}
};


/* Template implementation */

template <typename V>
void SplitRawSampleSet<V>::_read(fstream& file)
{
	const char delim = Base::Base::delim;
	const size_t nSkippedLines = Base::Base::nSkippedLines, headerLine = Base::Base::headerLine;
	
	// assume M makers and N samples
	// no headerLine
	
	marker::Set::ChromosomeMarkers* allMarkers = Base::Base::markers->unsortedChromosome();
	if (allMarkers == NULL) {
		throw logic_error("Markers should be populated on an unsorted chromosome prior to reading SplitRawSampleSet data");
	}
	
	marker::Set::ChromosomeMarkers::iterator markerIt = allMarkers->begin();
	marker::Set::ChromosomeMarkers::const_iterator markerEnd = allMarkers->end();
	
	// Use fileName without extension as sampleName
	string sampleName = name::filestem(Base::fileName);
	typename Base::RawSample* sample = Base::create(sampleName);
	
	size_t lineCount = 0;
	string line, markerName, chromName, discard;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines && lineCount != headerLine) {
			istringstream stream(line);
			
			// discard previous columns
			size_t colCount = 0;
			while (++colCount < dataColumn) {
				//stream >> discard;
				getline(stream, discard, delim);
			}
			
			if (markerIt == markerEnd) {
				throw runtime_error("Number of markers do not match the number of values for sample");
			}
			
			readSampleValue(stream, sample, markerIt->chromosome-1, delim);
			
			// next marker
			++markerIt;
		} else {
			// discard line
		}
	}
}

template <typename V> inline
void SplitRawSampleSet<V>::readSampleValue(istringstream& stream, typename Base::RawSample* sample, size_t chromIndex, const char delim) {
	typename Base::Value value;
	string s;
	getline(stream, s, delim);
	value = atof(s.c_str());
	sample->addToChromosome(chromIndex, value);
}


/* Template Specialization */

// Use the same specialization for AlleleSpecificCopyNumberValue and AlleleSpecificIntegerCopyNumberValue

#define SPECIALIZATION_TYPE AlleleSpecificCopyNumberValue
#include "SplitRawSampleSet.special"
#undef SPECIALIZATION_TYPE

#define SPECIALIZATION_TYPE AlleleSpecificIntegerCopyNumberValue
#include "SplitRawSampleSet.special"
#undef SPECIALIZATION_TYPE

#endif
