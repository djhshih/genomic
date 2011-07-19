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
	Segment(position startPos, position endPos, unsigned long numElements, float segValue) : start(startPos), end(endPos), nelem(numElements), value(segValue) {}
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


class SegmentedSampleSet;
class RawSampleSet;
class GenericSampleSet;

/*
template <typename V> class SegmentedSampleSet;
template <typename V> class RawSampleSet;
template <typename V> class GenericSampleSet;
*/

template <typename V>
class SampleSet
{
	friend class GenericSampleSet;
	//  for accessing the private clone() function
	//    and the read() and write() functions
public:
	typedef V Value;
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
class GenericSampleSet : public SampleSet<CopyNumberValue>
{
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

//template <typename V>
class RawSampleSet : public SampleSet<CopyNumberValue>
{
	friend class GenericSampleSet;
	friend class SegmentedSampleSet;
public:
	//typedef V Value;
	typedef LinearChromosome<Value> RawChromosome;
	typedef Sample<Value, RawChromosome > RawSample;
	typedef vector<RawSample*> Samples;
	typedef Samples::iterator SamplesIterator;
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
	RawSampleSet(marker::Set* markerSet) : SampleSet(markerSet) {}
	RawSampleSet(RawSampleSet& raw) {
		// TODO
	}
	RawSampleSet(SegmentedSampleSet& segmented);
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


//template <typename Value>
class SegmentedSampleSet : public SampleSet<CopyNumberValue>
{
	friend class GenericSampleSet;
	friend class RawSampleSet;
public:
	typedef LinearChromosome< Segment<Value> > SegmentedChromosome;
	typedef Sample< Segment<Value>, SegmentedChromosome > SegmentedSample;
	
	typedef vector<SegmentedSample*> Samples;
	typedef Samples::iterator SamplesIterator;
	typedef typename SegmentedSample::Chromosomes::iterator ChromosomesIterator;
	typedef SegmentedChromosome::iterator DataIterator;
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
	SegmentedSampleSet(marker::Set* markerSet) : SampleSet(markerSet) {}
	SegmentedSampleSet(SegmentedSampleSet& segmented) {
		byNames = segmented.byNames;
		// TODO
	}
	SegmentedSampleSet(RawSampleSet& raw);
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


class PicnicSampleSet : public RawSampleSet
{
	
};

class DchipSampleSet : public RawSampleSet
{
	
};

class CnagSampleSet : public RawSampleSet
{
	
};

#endif