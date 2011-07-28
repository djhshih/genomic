#ifndef genomic_Sample_h
#define genomic_Sample_h

#include "genomic.h"
#include "Tree.hpp"
#include "Marker.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>


//TODO use boost::spirit for parsing
//TODO allow use of input delimiter other than whitespace in all functions
//TODO improve sorting (enable skipping): used up 50GB of RAM on 3 picnic samples!


template <typename V>
class Segment
{
public:
	position start;
	position end;
	position nelem;
	V value;
	position length() {
		return end - start + 1;
	}
	Segment() {}
	Segment(position startPos, position endPos, unsigned long numElements, V segValue) : start(startPos), end(endPos), nelem(numElements), value(segValue) {}
	static bool compare(const Segment& a, const Segment& b) {
		return a.start < b.start; 
	}
};

template <typename T>
class Chromosome
{
public:
	typedef T DataType;
private:
	
public:
	size_t index;
	Chromosome() {}
	Chromosome(const size_t chromIndex) : index(chromIndex) {
	}
};

template <typename T>
class LinearChromosome : Chromosome<T>
{
public:
	typedef T DataType;
	typedef typename vector<T>::iterator iterator;
private:
	vector<T> items;
public:
	LinearChromosome() {}
	LinearChromosome(const size_t chromIndex) : Chromosome<T>(chromIndex) {
	}
	LinearChromosome(const LinearChromosome<T>& chr) 
	: Chromosome<T>(chr.index) {
		// deep copy
		for (size_t i = 0; i < chr.items.size(); ++i) {
			items.push_back(chr.items[i]);
		}
	}
	LinearChromosome<T>& operator=(const LinearChromosome<T>& chr) {
		LinearChromosome tmp(chr);
		swap(tmp);
		return *this;
	}
	LinearChromosome<T>& swap(LinearChromosome<T>& chr) {
		items.swap(chr.items);
	}
	~LinearChromosome() {
	}
	T& at(size_t i) {
		return items[i];
	}
	T& operator[](size_t i) {
		return items[i];
	}
	size_t size() {
		return items.size();
	}
	iterator begin() {
		return items.begin();
	}
	iterator end() {
		return items.end();
	}
	void push_back(const T& item) {
		items.push_back(item);
	}
	void pop_back() {
		items.pop_back();
	}
	void clear() {
		items.clear();
	}
};

template <typename T>
class KDTreeChromosome : Chromosome<T>
{
private:
	tree::KDTree< tree::Container< Segment<CopyNumberValue> > > items;
public:
	KDTreeChromosome() {
	}
};


template <typename Chromosome>
class Sample
{
public:
	typedef typename Chromosome::DataType T;
	typedef vector<Chromosome> Chromosomes;
	string name;
	string platform;
private:
	// vector of segments vectors for each Chromosome
	Chromosomes items;
public:
	Sample(const string& sampleName) : name(sampleName), items(nChromosomes) {
		// always allocate for all chromosomes
	}
	~Sample() {
		clear();
	}
	void clear() {
		items.clear();
	}
	Chromosome* operator[](size_t chromIndex) {
		if (chromIndex < items.size()) {
			return &(items[chromIndex]);
		}
		return NULL;
	}
	Chromosome* operator[](string& chromName) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) return &(items[k-1]);
		return NULL;
	}
	typename Chromosomes::iterator begin() {
		return items.begin();
	}
	typename Chromosomes::iterator end() {
		return items.end();
	}
	Chromosome* chromosome(size_t chromIndex) {
		if (chromIndex < items.size()) {
			return &(items[chromIndex]);
		}
		return NULL;
	}
	Chromosome* chromosome(const string& chromName) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) return &(items[k-1]);
		return NULL;
	}
	void addToChromosome(size_t chromIndex, const T& item) {
		if (chromIndex < items.size()) {
			items[chromIndex].push_back(item);
		}
	}
	void addToChromosome(const string& chromName, const T& item) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) {
			items[k-1].push_back(item);
		}
	}
	void remove(T item) {
		// not implemented
		// need to find segment in data structure
	}
	static bool compare(const Sample& a, const Sample& b) {
		return a.name < b.name;
	}
	static bool pcompare(Sample* a, Sample* b) {
		return a->name < b->name;
	}
};


template <typename V> class SegmentedSampleSet;
template <typename V> class RawSampleSet;
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
	void setIO(char _delim='\t', size_t _headerLine=1, size_t _nSkippedLines=0, bool _mergeSamples=false) {
		delim = _delim;
		headerLine = _headerLine;
		nSkippedLines = _nSkippedLines;
		mergeSamples = _mergeSamples;
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


// Generic sample set, chooses appropriately between possible types of sample set
// Use handle-body idiom
class GenericSampleSet : public SampleSet
{
public:
	typedef SampleSet Base;
private:
	// body representation
	// N.B. can only point to derived classes of SampleSet other than this class
	SampleSet* rep;
	
	SampleSet* clone() {
		return new GenericSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);
public:
	GenericSampleSet() : rep(NULL) {}
	GenericSampleSet(const GenericSampleSet& gset) {
		rep = gset.rep->clone();
	}
	~GenericSampleSet() {
		clear();
	}
	void clear() {
		delete rep;
		rep = NULL;
	}
	void sort() {
		rep->sort();
	}
	data::Type type() {
		return data::generic;
	}
};

template <typename V = CopyNumberValue>
class RawSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class SegmentedSampleSet<V>;
public:
	typedef SampleSet Base;
	typedef V Value;
	typedef LinearChromosome<Value> RawChromosome;
	typedef Sample<RawChromosome> RawSample;
	typedef vector<RawSample*> Samples;
	typedef typename Samples::iterator SamplesIterator;
	typedef typename RawSample::Chromosomes::iterator ChromosomesIterator;
	typedef typename RawChromosome::iterator DataIterator;

protected:
	Samples samples;

private:
	map<string, RawSample*> byNames;
	
	RawSampleSet* clone() {
		return new RawSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);

	void readSampleNames(istringstream& stream);
	void readSampleValues(istringstream& stream, size_t sampleStart, const string& chromName);

	void writeSampleNames(fstream& file, const char delim);
	void writeSampleValues(fstream& file, size_t chr, size_t markerIndex, const char delim);
	
public:
	RawSampleSet() {}
	RawSampleSet(marker::Set* markerSet) : SampleSet(markerSet) {}
	RawSampleSet(RawSampleSet& raw) {
		// TODO
	}
	RawSampleSet(SegmentedSampleSet<V>& segmented);
	~RawSampleSet() {
		clear();
	}
	data::Type type() {
		return data::raw;
	}
	void clear() {
		SamplesIterator i, end = samples.end();
		for (i = samples.begin(); i != end; ++i) {
			// delete object pointed to by Sample* pointer
			delete (*i);
		}
		samples.clear();
		byNames.clear();
		marker::manager.unref(markers);
	}
	RawSample* create(const string& sampleName) {
		RawSample* sam = byNames[sampleName];
		if (sam == NULL) {
			// sample does not exist: create it
			sam = new RawSample(sampleName);
			samples.push_back(sam);
			// register name
			byNames[sampleName] = sam;
		}
		return sam;
	}
	void sort();
};


template <typename V = CopyNumberValue>
class SegmentedSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class RawSampleSet<V>;
public:
	typedef SampleSet Base;
	typedef V Value;
	typedef LinearChromosome< Segment<Value> > SegmentedChromosome;
	typedef Sample<SegmentedChromosome> SegmentedSample;
	
	typedef vector<SegmentedSample*> Samples;
	typedef typename Samples::iterator SamplesIterator;
	typedef typename SegmentedSample::Chromosomes::iterator ChromosomesIterator;
	typedef typename SegmentedChromosome::iterator DataIterator;
private:
	Samples samples;
	map<string, SegmentedSample*> byNames;
	
	SegmentedSampleSet* clone() {
		return new SegmentedSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);
	
	void readSegment(fstream& file, Segment<V>& seg) {
		file >> seg.start >> seg.end >> seg.nelem >> seg.value;
	}
	
	size_t _find(SegmentedChromosome& array, position x);
	
public:
	SegmentedSampleSet() {}
	SegmentedSampleSet(marker::Set* markerSet) : SampleSet(markerSet) {}
	SegmentedSampleSet(SegmentedSampleSet& segmented) {
		byNames = segmented.byNames;
		// TODO
	}
	SegmentedSampleSet(RawSampleSet<V>& raw);
	~SegmentedSampleSet() {
		clear();
	}
	data::Type type() {
		return data::segmented;
	}
	void clear() {
		SamplesIterator i, end = samples.end();
		for (i = samples.begin(); i != end; ++i) {
			// delete object pointed to by Sample* pointer
			delete (*i);
		}
		samples.clear();
		byNames.clear();
	}
	SegmentedSample* create(const string& sampleName) {
		SegmentedSample* sam = byNames[sampleName];
		if (sam == NULL) {
			// sample does not exist: create it
			sam = new SegmentedSample(sampleName);
			samples.push_back(sam);
			// register name
			byNames[sampleName] = sam;
		}
		return sam;
	}
	void sort();
	
	size_t find(const string& sampleName, size_t chromIndex, position start) {
		return _find(*(byNames[sampleName]->chromosome(chromIndex)), start);
	}
	size_t find(SegmentedSample* sample, size_t chromIndex, position start) {
		return _find(*(sample->chromosome(chromIndex)), start);
	}
	
	void filter(SegmentedSampleSet& ref, float diceThreshold);
};

template <typename V>
class ReferenceSegmentedSampleSet : public SegmentedSampleSet<V>
{
public:
	typedef SegmentedSampleSet<V> Base;
public:
	ReferenceSegmentedSampleSet() {
		Base::Base::mergeSamples = true;
	}
};

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

class PicnicSampleSet : public SplitRawSampleSet<AlleleSpecificIntegerCopyNumberValue>
{
public:
	typedef SplitRawSampleSet<AlleleSpecificIntegerCopyNumberValue> Base;
public:
	// set dataColumn to 5
	PicnicSampleSet() : Base(5) {
		Base::Base::delim = ',';
		Base::Base::headerLine = 0;
	}
	data::Type type() {
		return data::picnic;
	}
};

class DchipSampleSet : public RawSampleSet<CopyNumberValue>
{
	
};

class CnagSampleSet : public RawSampleSet<CopyNumberValue>
{
	
};



/* Template implementation */
/* Required to be in the same file as the definitions */

template <typename V>
RawSampleSet<V>::RawSampleSet(SegmentedSampleSet<V>& set)
{
	//TODO
	// N.B. throw error if set.markers does not exist
	//   Marker information in set is required for creating RawSampleSet
}

template <typename V>
void RawSampleSet<V>::_read(fstream& file)
{
	const char delim = Base::delim;
	const size_t nSkippedLines = Base::nSkippedLines, headerLine = Base::headerLine;
	marker::Set* markers = Base::markers;
	
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	
	bool readMarkers = (markers->empty()) ? true : false;
	
	size_t lineCount = 0;
	size_t sampleStart = samples.size()-1; 
	string line, markerName, chromName, discard;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines) {
			if (lineCount == headerLine) {
				istringstream stream(line);
				// discard the marker information columns (3)
				stream >> discard >> discard >> discard;
				
				readSampleNames(stream);
			} else {
				istringstream stream(line);
				position pos;
				stream >> markerName >> chromName >> pos;
				
				if (readMarkers) {
					size_t chr = mapping::chromosome[chromName];
					// ignore unknown chromosome: continue to next line
					if (chr == 0) continue;
					// create marker
					marker::Marker marker(markerName, chr, pos);
					markers->at(chr-1).push_back(marker);
				}
				
				readSampleValues(stream, sampleStart, chromName);
			}
		} else {
			// discard line
		}
	}
}

template <typename V> inline
void RawSampleSet<V>::readSampleNames(istringstream& stream) {
	while (!stream.eof()) {
		string sampleName;
		stream >> sampleName;
		// create sample with $sampleName	
		create(sampleName);
	}
}

template <typename V> inline
void RawSampleSet<V>::readSampleValues(istringstream& stream, size_t sampleStart, const string& chromName) {
	size_t i = sampleStart;
	Value value;
	while (!stream.eof()) {
		stream >> value;
		// create point at specified chromosome
		samples[++i]->addToChromosome(chromName, value);
	}
}

template <typename V>
void RawSampleSet<V>::_write(fstream& file)
{
	const char delim = Base::delim;
	marker::Set* markers = Base::markers;
	
	file << "marker" << delim << "chromosome" << delim << "position";
	
	writeSampleNames(file, delim);
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers->size(); ++chr) {
		for (size_t markerIndex = 0; markerIndex < markers->at(chr).size(); ++markerIndex) {
			
			// print marker information
			file << markers->at(chr)[markerIndex].name << delim << markers->at(chr)[markerIndex].chromosome << delim << markers->at(chr)[markerIndex].pos;
			
			writeSampleValues(file, chr, markerIndex, delim);
		}
	}
}

template <typename V> inline
void RawSampleSet<V>::writeSampleNames(fstream& file, const char delim) {
	// print sample names
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << endl;
}

template <typename V> inline
void RawSampleSet<V>::writeSampleValues(fstream& file, size_t chr, size_t markerIndex, const char delim) {
	// iterate through samples to print values, selected the specified chromosome and marker
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it)[chr]->at(markerIndex);
	}
	file << endl;
}

template <typename V>
void RawSampleSet<V>::sort()
{
	marker::Set* markers = Base::markers;
	
	// Construct order vector for obtaining a sorted index of markers
	vector< pair<position, size_t> > order;
	
	for (size_t chri = 0; chri < markers->size(); ++chri) {	
		
		// Construct order vector for obtaining a sorted index of markers
		// Additionally, replicate the markers on the curent chromosome;
		//               replicate the chromosome for all samples
		marker::Set::ChromosomeMarkers chromosomeMarkers;
		vector<RawChromosome> samplesChromosomeCopy;
		vector< pair<position, size_t> > order;
		for (size_t j = 0; j < markers->at(chri).size(); ++j) {
			order.push_back( make_pair(markers->at(chri)[j].pos, j) );
			
			chromosomeMarkers.push_back(markers->at(chri)[j]);
			
			for (size_t s = 0; s < samples.size(); ++s) {
				samplesChromosomeCopy.push_back(*(samples[s]->chromosome(chri)));
			}
		}
		
		// Sort on the order vector instead of the original vector<Markers*>,
		//   in order to obtain the sorted index
		std::sort(order.begin(), order.end(), &compare::pair<position, size_t>);
		
		// Now, order is sorted by genomic position, and order[i].second
		//   contains each sorted index
		
		// Sort the markers on current chromosome, and each sample
		for (size_t j = 0; j < markers->at(chri).size(); ++j) {
			size_t index = order[j].second;
			// Set the marker to the corresponding sorted marker
			markers->at(chri)[j] = chromosomeMarkers[index];
			
			// Iterate through samples, set the corresponding marker in the current chromosome
			for (size_t s = 0; s < samples.size(); ++s) {
				samples[s]->chromosome(chri)->at(j) = samplesChromosomeCopy[s][index];
			}
		}
		
	}
}

template <typename V>
SegmentedSampleSet<V>::SegmentedSampleSet(RawSampleSet<V>& raw)
{
	clear();
	// use iterators to avoid assuming RawSampleSet stores data in vectors
	// however, need to assume that markers are stored in vectors, for looking up marker information
	
	typedef typename RawSampleSet<V>::Samples::iterator RawSamplesIterator;
	typedef typename RawSampleSet<V>::ChromosomesIterator RawChromosomesIterator;
	typedef typename RawSampleSet<V>::DataIterator RawDataIterator;
	
	// iterate through samples in $raw
	RawSamplesIterator it, end = raw.samples.end();
	for (it = raw.samples.begin(); it != end; ++it) {
		// create sample
		SegmentedSample* sam = create((*it)->name);
		// iterate through chromosome in sample
		size_t chr = 0;
		RawChromosomesIterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through markers on chromosome
			// CANNOT assume markers are stored in a vector
			RawDataIterator markerIt = chrIt->begin(), markerEnd = chrIt->end();
			position markerIndex, startMarkerIndex = 0;
			V prevValue;
			//Value prevValue;
			if (markerIt != markerEnd) {
				prevValue = *markerIt;
				// start at second marker
				++markerIt;
				markerIndex = 1;
				while (markerIt != markerEnd) {
					if (!eq(*markerIt, prevValue)) {
						// segment ended: store segment from $startMarkerIndex to $markerIndex-1
						Segment<Value> seg(
							raw.markers->at(chr)[startMarkerIndex].pos,
							raw.markers->at(chr)[markerIndex-1].pos,
							(markerIndex-1) - startMarkerIndex + 1,
							prevValue
						);
						sam->chromosome(chr)->push_back(seg);
						// start new segment
						startMarkerIndex = markerIndex;
						prevValue = *markerIt;
					}
					++markerIndex;
					++markerIt;
				}
				// store last segment
				// handling is same whether last segment is the last marker alone or
				// 	laste segment ends on the last marker
				Segment<Value> seg(
					raw.markers->at(chr)[startMarkerIndex].pos,
					raw.markers->at(chr)[markerIndex-1].pos,
					(markerIndex-1) - startMarkerIndex + 1,
					prevValue
				);
				sam->chromosome(chr)->push_back(seg);
			}
			++chr;
		}  // for chrIt
	} // for it
}

template <typename V>
void SegmentedSampleSet<V>::_read(fstream& file)
{
	const char delim = Base::delim;
	const size_t nSkippedLines = Base::nSkippedLines, headerLine = Base::headerLine;
	
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	
	size_t lineCount = 0;
	string line, sampleName, chromName;
	Segment<Value>* seg;
	while (true) {
		if (++lineCount > nSkippedLines && lineCount != headerLine) {
			file >> sampleName >> chromName;
			if (file.eof()) break;
			// ignore unknown chromosome: continue to next line
			if (mapping::chromosome[chromName] == 0) continue;
			// create segment at specified chromosome
			Segment<V> seg;
			readSegment(file, seg);
			create(sampleName)->addToChromosome(chromName, seg);
			//trace("%s %s %d %d %d %f\n", sampleName.c_str(), chromName.c_str(), seg.start, seg.end, seg.nelem, seg.value);
		} else {
			// discard line
			getline(file, line);
		}
	}
}

template <typename V>
void SegmentedSampleSet<V>::_write(fstream& file)
{
	const char delim = Base::delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "state" << endl;
	
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->nelem << delim << segIt->value << endl;
			}
			++chr;
		}
	}
}

template <typename V>
void SegmentedSampleSet<V>::sort()
{
	// Sort samples by name
	std::sort(samples.begin(), samples.end(), &SegmentedSample::pcompare);
	// Iterate through samples and chromosomes therefore, sort segments
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			std::sort(chrIt->begin(), chrIt->end(), &Segment<Value>::compare);
		}
	}
}

template <typename V>
size_t SegmentedSampleSet<V>::_find(SegmentedChromosome& array, position x) {
	// Initialize left and right beyond array bounds
	size_t left = -1, right = array.size();
	while (left + 1 != right) {
		// Check middle of remaining subarray
		size_t i = (left + right) / 2;
		if (x < array[i].start) right = i;      // in the left half
		if (x == array[i].start) return i;      // found
		if (x > array[i].start) left = i;       // in the left half
	}
	// x is not found in the array
	// return index of the value that is greatest value lower than the query
	return ( (left == -1) ? 0 : left );
}

template <typename V>
void SegmentedSampleSet<V>::filter(SegmentedSampleSet& ref, float diceThreshold)
{
	Samples oldSamples;
	// iterate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		ChromosomesIterator chrIt;
		const ChromosomesIterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			DataIterator segIt;
			const DataIterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				// Compare against segments on reference
				SamplesIterator refIt;
				const SamplesIterator refEnd = ref.samples.end();
				bool filterSegment = false;
				for (refIt = ref.samples.begin(); refIt != refEnd; ++refIt) {
					// determine lower and upper bounds
					SegmentedChromosome* refChrom = (**refIt)[chri];
					position_diff lower = 2*(diceThreshold-1)/diceThreshold*segIt->end + (2-diceThreshold)/diceThreshold*(segIt->start - 1) + 1;
					if (lower < 0) lower = 0;
					position_diff upper = 2*(1-diceThreshold)/(2-diceThreshold)*segIt->end + diceThreshold/(2-diceThreshold)*(segIt->start - 1) + 1;
					//cout << segIt->start << " " << segIt->end << " " << lower << " " << upper << endl;
					size_t lowerIndex = ref.find(*refIt, chri, lower);
					size_t upperIndex = ref.find(*refIt, chri, upper) + 1;
					if (upperIndex >= refChrom->size()) upperIndex = refChrom->size()-1;
					//size_t lowerIndex = 0, upperIndex = refChrom->size()-1;
					//cout << "Index: " << lowerIndex << ", " << upperIndex << endl;
					for (size_t i = lowerIndex; i <= upperIndex; ++i) {
						// calculate Dice coefficient
						position_diff intersection = min(refChrom->at(i).end, segIt->end) - max(refChrom->at(i).start, segIt->start) + 1;
						//cout << segIt->start << " " << refChrom->at(i).start << " " << intersection << endl;
						if (intersection > 0) {
							float dice = 2 * float(intersection) / (refChrom->at(i).length() + segIt->length());
							if (dice > diceThreshold) {
								//cout << "Filter: " << segIt->start << " " << refChrom->at(i).start << " " << dice << endl;
								// Mark segment for deletion
								segIt->nelem = 0;
								filterSegment = true;
								break;
							}
						}
					}
					if (filterSegment) break;
				}
			}
			++chri;
		}
		
		// copy current samples
		oldSamples.push_back(*it);
	}
	
	// Create new sample set with marked segments removed
	samples.clear();
	byNames.clear();
	end = oldSamples.end();
	for (it = oldSamples.begin(); it != end; ++it) {
		SegmentedSample* sample = create((*it)->name);
		ChromosomesIterator chrIt;
		const ChromosomesIterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			DataIterator segIt;
			const DataIterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				if (segIt->nelem > 0) {
					Segment<V> seg(segIt->start, segIt->end, segIt->nelem, segIt->value);
					sample->addToChromosome(chri, seg);
				}
			}
			++chri;
		}
		// delete old sample
		delete (*it);
	}
	oldSamples.clear();
	
}


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
	
	marker::Set::MarkersIterator markerIt = allMarkers->begin();
	const marker::Set::MarkersIterator markerEnd = allMarkers->end();
	
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


//// Template Specialization

// Use the same specialization for AlleleSpecificCopyNumberValue and AlleleSpecificIntegerCopyNumberValue

#define SPECIALIZATION_TYPE AlleleSpecificCopyNumberValue
#include "Sample_AlleleSpecific_specialization.h"
#undef SPECIALIZATION_TYPE

#define SPECIALIZATION_TYPE AlleleSpecificIntegerCopyNumberValue
#include "Sample_AlleleSpecific_specialization.h"
#undef SPECIALIZATION_TYPE


#endif
