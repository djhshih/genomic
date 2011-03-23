#ifndef genomic_Sample_h
#define genomic_Sample_h

#include "genomic.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

class Segment
{
public:
	position start;
	position end;
	position nelem;
	float value;
	position length() {
		return end - start + 1;
	}
	Segment() {}
	Segment(position startPos, position endPos, unsigned long numElements, float segValue) : start(startPos), end(endPos), nelem(numElements), value(segValue) {}
};

template<typename T>
class Sample
{
	friend class SampleSet;
	friend class SegmentedSampleSet;
	friend class RawSampleSet;
public:
	typedef vector<T*> Chromosome;
	typedef vector<Chromosome> Chromosomes;
	string name;
	string platform;
private:
	// vector of segments vectors for each Chromosome
	Chromosomes items;
public:
	Sample(const string& sampleName) : name(sampleName) {
		items.resize(nChromosomes);
	}
	~Sample() {
		for (typename Chromosomes::iterator i = items.begin(); i != items.end(); ++i) {
			for (typename Chromosome::iterator j = i->begin(); j != i->end(); ++j) {
				// delete object pointed to by Segment* pointer
				delete (*j);
			}
		}
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
	T* createAt(size_t chromIndex) {
		if (chromIndex < items.size()) {
			T* item = new T();
			items[chromIndex].push_back(item);
			return item;
		}
		return NULL;
	}
	T* createAt(string& chromName) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) {
			T* item = new T();
			items[k-1].push_back(item);
			return item;
		}
		return NULL;
	}
	void remove(T* item) {
		// not implemented
		// need to find segment in data structure
	}
};


class SegmentedSampleSet;
class RawSampleSet;

class SampleSet
{
public:
	SampleSet() {}
	//virtual SampleSet(SampleSet& set) = 0;
	virtual ~SampleSet() {
		if (file.is_open()) file.close();
	}
	virtual void read(const string& fileName) = 0;
	virtual void write(const string& fileName) = 0;
	virtual void clear() = 0;
	virtual data::Type type() = 0;
protected:
	//virtual SampleSet* cloneFrom(SampleSet* set) = 0;
	static const char delim;
	fstream file;
};


// Generic sample set, chooses appropriately between possible types of sample set
// Use handle-body idiom
class GenericSampleSet : public SampleSet
{
private:
	// body representation
	// cannot only point to derived classes of SampleSet other than this class
	SampleSet* rep;
	/*SampleSet* cloneFrom(SampleSet* set) {
		return rep->cloneFrom(set);
	}*/
public:
	GenericSampleSet() : rep(NULL) {}
	~GenericSampleSet() {
		clear();
	}
	void read(const string& fileName);
	void write(const string& fileName);
	void clear() {
		delete rep;
		rep = NULL;
	}
	data::Type type() {
		return data::generic;
	}
};


//TODO Ensure SegmentedSampleSet and RawSampleSet sort input after storing


class RawSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class SegmentedSampleSet;
public:
	typedef float Value;
	typedef vector<Sample<Value>*> Samples;
	typedef Samples::iterator SamplesIterator;
	typedef Sample<Value>::Chromosomes::iterator ChromosomesIterator;
	typedef Sample<Value>::Chromosome::iterator DataIterator;
	// vector of vector of markers organized by chromosomes
	typedef vector< vector<Marker*> > Markers;
private:
	Samples samples;
	Markers markers;
	map<string, Sample<Value>*> byNames;
	
	/*RawSampleSet* cloneFrom(SampleSet* set) {
		return new RawSampleSet(*set);
	}*/
public:
	RawSampleSet() {
		markers.resize(nChromosomes);
	}
	RawSampleSet(SampleSet& set);
	RawSampleSet(RawSampleSet& raw);
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
		Markers::iterator mi, mend = markers.end();
		for (mi = markers.begin(); mi != mend; ++mi) {
			vector<Marker*>::iterator mj, mjend = mi->end();
			for (mj = mi->begin(); mj != mjend; ++mj) {
				delete (*mj);
			}
		}
		markers.clear();
		markers.resize(nChromosomes);
	}
	void read(const string& fileName);
	void write(const string& fileName);
	Sample<Value>* sample(const string& sampleName) {
		Sample<Value>* sam = byNames[sampleName];
		if (sam == NULL) {
			// sample does not exist: create it
			sam = new Sample<Value>(sampleName);
			samples.push_back(sam);
			// register name
			byNames[sampleName] = sam;
		}
		return sam;
	}
};



class SegmentedSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class RawSampleSet;
public:
	typedef vector<Sample<Segment>*> Samples;
	typedef Samples::iterator SamplesIterator;
	typedef Sample<Segment>::Chromosomes::iterator ChromosomesIterator;
	typedef Sample<Segment>::Chromosome::iterator DataIterator;
private:
	Samples samples;
	map<string, Sample<Segment>*> byNames;
	
	/*SegmentedSampleSet* cloneFrom(SampleSet* set) {
		return new SegmentedSampleSet(*set);
	}*/
public:
	SegmentedSampleSet() {
	}
	SegmentedSampleSet(SegmentedSampleSet& segmented);
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
	void read(const string& fileName);
	void write(const string& fileName);
	Sample<Segment>* sample(const string& sampleName) {
		Sample<Segment>* sam = byNames[sampleName];
		if (sam == NULL) {
			// sample does not exist: create it
			sam = new Sample<Segment>(sampleName);
			samples.push_back(sam);
			// register name
			byNames[sampleName] = sam;
		}
		return sam;
	}
};




#endif