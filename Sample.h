#ifndef genomic_Sample_h
#define genomic_Sample_h

#include "genomic.h"
#include "Tree.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

typedef float CopyNumberValue;

struct AlleleSpecificCopyNumberValue
{
	float a, b;
	bool operator!=(const AlleleSpecificCopyNumberValue& y) {
		return !(a == y.a && b == y.b);
	}
	bool operator==(const AlleleSpecificCopyNumberValue& y) {
		return a == y.a && b == y.b;
	}
};

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
		// could not use iterator here... ambiguity issue with T?
		/*
		iterator i, end = chr.items.end();
		for (i = chr.items.begin(); i != end; ++i) {
			items.push_back(*i);
		}
		*/
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

// TODO let Sample determine type from Chromosome

template <typename T, typename Chromosome>
class Sample
{
public:
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
	Chromosome* chromosome(string& chromName) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) return &(items[k-1]);
		return NULL;
	}
	void addToChromosome(size_t chromIndex, T& item) {
		if (chromIndex < items.size()) {
			items[chromIndex].push_back(item);
		}
	}
	void addToChromosome(string& chromName, T& item) {
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

/*
class SegmentedSampleSet;
class RawSampleSet;
class GenericSampleSet;
*/

template <typename V> class SegmentedSampleSet;
template <typename V> class RawSampleSet;
template <typename V> class GenericSampleSet;

//template <typename V>
template <typename V = CopyNumberValue>
class SampleSet
{
	friend class GenericSampleSet<V>;
	//  for accessing the private clone() function
	//    and the read() and write() functions
public:
	//typedef V Value;
public:
	SampleSet() : markers(NULL) {}
	SampleSet(marker::Set* markerSet) : markers(markerSet) {}
	//virtual SampleSet(SampleSet& set) = 0;
	virtual ~SampleSet() {
		if (file.is_open()) file.close();
	}
	void read(const string& fileName) {
		// use fileName also as platform name
		read(fileName, fileName);
	}
	void read(const string& fileName, const string& platform) {
		file.open(fileName.c_str(), ios::in);
		if (!file.is_open()) throw runtime_error("Failed to open input file");
		read(file, platform, fileName);
		file.close();
		trace("Read file %s\n", fileName.c_str());
	}
	void read(fstream& file, const string& platform, const string& fileType) {
		clear();
		markers = marker::manager.create(platform);
		this->fileName = "." + fileType;
		
		_read(file);
		sort();
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
	static const char delim;
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
//template <typename V>
template <typename V = CopyNumberValue>
class GenericSampleSet : public SampleSet<V>
{
public:
	typedef SampleSet<V> Base;
private:
	// body representation
	// N.B. can only point to derived classes of SampleSet other than this class
	SampleSet<V>* rep;
	
	SampleSet<V>* clone() {
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

//template <typename V>
template <typename V = CopyNumberValue>
class RawSampleSet : public SampleSet<V>
{
	friend class GenericSampleSet<V>;
	friend class SegmentedSampleSet<V>;
public:
	typedef SampleSet<V> Base;
	typedef V Value;
	typedef LinearChromosome<Value> RawChromosome;
	typedef Sample<Value, RawChromosome > RawSample;
	typedef vector<RawSample*> Samples;
	typedef typename Samples::iterator SamplesIterator;
	typedef typename RawSample::Chromosomes::iterator ChromosomesIterator;
	typedef typename RawChromosome::iterator DataIterator;
private:
	Samples samples;
	map<string, RawSample*> byNames;
	
	RawSampleSet* clone() {
		return new RawSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);
public:
	RawSampleSet() {}
	RawSampleSet(marker::Set* markerSet) : SampleSet<V>(markerSet) {}
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
		marker::manager.unref(SampleSet<V>::markers);
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


//template <typename V = CopyNumberValue>
template <typename V>
class SegmentedSampleSet : public SampleSet<V>
{
	friend class GenericSampleSet<V>;
	friend class RawSampleSet<V>;
public:
	typedef SampleSet<V> Base;
	typedef V Value;
	typedef LinearChromosome< Segment<Value> > SegmentedChromosome;
	typedef Sample< Segment<Value>, SegmentedChromosome > SegmentedSample;
	
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
public:
	SegmentedSampleSet() {}
	SegmentedSampleSet(marker::Set* markerSet) : SampleSet<V>(markerSet) {}
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
	
	void filter(SegmentedSampleSet& ref);
};


class PicnicSampleSet : public RawSampleSet<AlleleSpecificCopyNumberValue>
{
	
};

class DchipSampleSet : public RawSampleSet<CopyNumberValue>
{
	
};

class CnagSampleSet : public RawSampleSet<CopyNumberValue>
{
	
};



/* Template implementation */
/* Required to be in the same file as the definitions */

template <typename V> const char SampleSet<V>::delim = '\t';

template <typename V>
void GenericSampleSet<V>::_read(fstream& file)
{
	const string& fileName = Base::fileName;
	marker::Set* markers = Base::markers;
	
	string ext = fileName.substr(fileName.find_last_of('.')+1);
	switch (mapping::extension[ext]) {
		case data::raw:
			rep = new RawSampleSet<V>(markers);
			break;
		case data::segmented:
			rep = new SegmentedSampleSet<V>(markers);
			break;
		default:
			throw runtime_error("Cannot determine file type from file name extension");
	}
	rep->_read(file);
}

template <typename V>
void GenericSampleSet<V>::_write(fstream& file)
{
	const string& fileName = Base::fileName;
	
	string ext = fileName.substr(fileName.find_last_of('.')+1);
	
	// cast $rep to appropriate type
	// runtime checking is skipped (i.e. use static_cast instead of dynamic_case)
	//  since exact type can be determined
	SampleSet<V>* tmp = NULL;
	switch (mapping::extension[ext]) {
		case data::raw:
			switch (rep->type()) {
				case data::raw:
					// nothing needs to be done
					break;
				case data::segmented:
					tmp = new RawSampleSet<V>(*static_cast<SegmentedSampleSet<V>*>(rep));
					swap(rep, tmp);
					delete tmp;
					break;
			}
			break;
		case data::segmented:
			switch (rep->type()) {
				case data::raw:
					tmp = new SegmentedSampleSet<V>(*static_cast<RawSampleSet<V>*>(rep));
					swap(rep, tmp);
					delete tmp;
					break;
				case data::segmented:
					// nothing to be done
					break;
			}
			break;
		default:
			throw runtime_error("Cannot determine file type from file name extension");
	}
	rep->_write(file);
}

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
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	marker::Set* markers = Base::markers;
	
	string line;
	
	size_t nSkippedLines = 0;
	size_t headerLine = 1;
	size_t lineCount = 0;
	string markerName, chromName, sampleName, discard;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines) {
			if (lineCount == headerLine) {
				istringstream stream(line);
				// discard the marker information columns (3)
				stream >> discard >> discard >> discard;
				// create samples
				while (!stream.eof()) {
					stream >> sampleName;
					// create sample with $sampleName	
					create(sampleName);
				}
			} else {
				istringstream stream(line);
				position pos;
				stream >> markerName >> chromName >> pos;
				size_t chr = mapping::chromosome[chromName];
				// ignore unknown chromosome: continue to next line
				if (chr == 0) continue;
				// create marker
				marker::Marker marker(markerName, chr, pos);
				markers->at(chr-1).push_back(marker);
				
				Value value;
				size_t i = -1;
				while (!stream.eof()) {
					stream >> value;
					// create point at specified chromosome
					samples[++i]->addToChromosome(chromName, value);
				}
			}
		} else {
			// discard line
		}
	}
}

template <typename V>
void RawSampleSet<V>::_write(fstream& file)
{
	const char delim = Base::delim;
	marker::Set* markers = Base::markers;
	
	file << "marker" << delim << "chromosome" << delim << "position";
	
	// print sample names
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << endl;
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers->size(); ++chr) {
		for (unsigned long markerIndex = 0; markerIndex < markers->at(chr).size(); ++markerIndex) {
			
			// print marker information
			file << markers->at(chr)[markerIndex].name << delim << markers->at(chr)[markerIndex].chromosome << delim << markers->at(chr)[markerIndex].pos;
			
			// iterate through samples to print values, selected the specified chromosome and marker
			SamplesIterator it, end = samples.end();
			for (it = samples.begin(); it != end; ++it) {
				file << delim << (**it)[chr]->at(markerIndex);
			}
			file << endl;
			
		}
	}
}

template <> inline
void RawSampleSet<AlleleSpecificCopyNumberValue>::_read(fstream& file)
{
	/*
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	marker::Set* markers = Base::markers;
	
	string line;
	
	size_t nSkippedLines = 0;
	size_t headerLine = 1;
	size_t lineCount = 0;
	string markerName, chromName, sampleName, discard;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines) {
			if (lineCount == headerLine) {
				istringstream stream(line);
				// discard the marker information columns (3)
				stream >> discard >> discard >> discard;
				// create samples
				while (!stream.eof()) {
					stream >> sampleName;
					// create sample with $sampleName	
					create(sampleName);
				}
			} else {
				istringstream stream(line);
				position pos;
				stream >> markerName >> chromName >> pos;
				size_t chr = mapping::chromosome[chromName];
				// ignore unknown chromosome: continue to next line
				if (chr == 0) continue;
				// create marker
				marker::Marker marker(markerName, chr, pos);
				markers->at(chr-1).push_back(marker);
				
				Value value;
				size_t i = -1;
				while (!stream.eof()) {
					stream >> value;
					// create point at specified chromosome
					samples[++i]->addToChromosome(chromName, value);
				}
			}
		} else {
			// discard line
		}
	}
	*/
}

template <> inline
void RawSampleSet<AlleleSpecificCopyNumberValue>::_write(fstream& file)
{
	/*
	const char delim = Base::delim;
	marker::Set* markers = Base::markers;
	
	file << "marker" << delim << "chromosome" << delim << "position";
	
	// print sample names
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << endl;
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers->size(); ++chr) {
		for (unsigned long markerIndex = 0; markerIndex < markers->at(chr).size(); ++markerIndex) {
			
			// print marker information
			file << markers->at(chr)[markerIndex].name << delim << markers->at(chr)[markerIndex].chromosome << delim << markers->at(chr)[markerIndex].pos;
			
			// iterate through samples to print values, selected the specified chromosome and marker
			SamplesIterator it, end = samples.end();
			for (it = samples.begin(); it != end; ++it) {
				file << delim << (**it)[chr]->at(markerIndex);
			}
			file << endl;
			
		}
	}
	*/
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
					if (*markerIt != prevValue) {
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
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	
	string line;
	size_t nSkippedLines = 1;
	size_t lineCount = 0;
	string sampleName, chromName;
	Segment<Value>* seg;
	while (true) {
		if (++lineCount > nSkippedLines) {
			file >> sampleName >> chromName;
			if (file.eof()) break;
			// ignore unknown chromosome: continue to next line
			if (mapping::chromosome[chromName] == 0) continue;
			// create segment at specified chromosome
			Segment<Value> seg;
			file >> seg.start >> seg.end >> seg.nelem >> seg.value;
			create(sampleName)->addToChromosome(chromName, seg);
			//trace("%s %s %d %d %d %f\n", sampleName.c_str(), chromName.c_str(), seg->start, seg->end, seg->nelem, seg->value);
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
	
	// iteratrate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		// can assume that Sample are stored chromosomes in a vector
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			// CANNOT assume segments are stored in a vector
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->nelem << delim << segIt->value << endl;
			}
			++chr;
		}
	}
}

template <> inline
void SegmentedSampleSet<AlleleSpecificCopyNumberValue>::_read(fstream& file)
{
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	
	string line;
	size_t nSkippedLines = 1;
	size_t lineCount = 0;
	string sampleName, chromName;
	Segment<Value>* seg;
	while (true) {
		if (++lineCount > nSkippedLines) {
			file >> sampleName >> chromName;
			if (file.eof()) break;
			// ignore unknown chromosome: continue to next line
			if (mapping::chromosome[chromName] == 0) continue;
			// create segment at specified chromosome
			Segment<Value> seg;
			file >> seg.start >> seg.end >> seg.nelem >> seg.value.a >> seg.value.b;
			create(sampleName)->addToChromosome(chromName, seg);
		} else {
			// discard line
			getline(file, line);
		}
	}
}

template <> inline
void SegmentedSampleSet<AlleleSpecificCopyNumberValue>::_write(fstream& file)
{
	const char delim = Base::delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "stateA" << delim << "stateB" << endl;
	
	// iteratrate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		// can assume that Sample are stored chromosomes in a vector
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			// CANNOT assume segments are stored in a vector
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->nelem << delim << segIt->value.a << delim << segIt->value.b << endl;
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
void SegmentedSampleSet<V>::filter(SegmentedSampleSet& ref)
{
	// iteratrate through samples
	SamplesIterator it;
	const SamplesIterator end = samples.end();
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
				//TODO
			}
			++chri;
		}
	}
}


#endif