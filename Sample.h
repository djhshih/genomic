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
	static bool compare(const Segment& a, const Segment& b) {
		return a.start < b.start; 
	}
};

template<typename T>
class Chromosome
{
private:
	
public:
	size_t index;
	Chromosome() {}
	Chromosome(const size_t chromIndex) : index(chromIndex) {
	}
};

template<typename T>
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

template<typename T>
class KDTreeChromosome : Chromosome<T>
{
private:
	tree::KDTree< tree::Container<Segment> > items;
public:
	KDTreeChromosome() {
	}
};

template<typename T, typename Chromosome>
class Sample
{
public:
	//typedef vector<T*> Chromosome;
	typedef vector<Chromosome> Chromosomes;
	string name;
	string platform;
private:
	// vector of segments vectors for each Chromosome
	Chromosomes items;
public:
	Sample(const string& sampleName) : name(sampleName) {
		// always allocate for all chromosomes
		items.resize(nChromosomes);
	}
	~Sample() {
		clear();
	}
	void clear() {
		items.clear();
		/*
		for (typename Chromosomes::iterator i = items.begin(); i != items.end(); ++i) {
			for (typename Chromosome::iterator j = i->begin(); j != i->end(); ++j) {
				// delete object pointed to by Segment* pointer
				delete (*j);
			}
		}
		*/
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

class SampleSet
{
	friend class GenericSampleSet;
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
	static const char delim;
	fstream file;
private:
	virtual SampleSet* clone() = 0;
};


// Generic sample set, chooses appropriately between possible types of sample set
// Use handle-body idiom
class GenericSampleSet : public SampleSet
{
private:
	// body representation
	// N.B. can only point to derived classes of SampleSet other than this class
	SampleSet* rep;
	
	SampleSet* clone() {
		return new GenericSampleSet(*this);
	}
public:
	GenericSampleSet() : rep(NULL) {}
	GenericSampleSet(const GenericSampleSet& gset) {
		rep = gset.rep->clone();
	}
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


class RawSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class SegmentedSampleSet;
public:
	typedef float Value;
	typedef LinearChromosome<Value> RawChromosome;
	typedef Sample<Value, RawChromosome > RawSample;
	typedef vector<RawSample*> Samples;
	typedef Samples::iterator SamplesIterator;
	typedef typename RawSample::Chromosomes::iterator ChromosomesIterator;
	typedef typename RawChromosome::iterator DataIterator;
	// vector of vector of markers organized by chromosomes
	typedef vector<Marker*> ChromosomeMarkers;
	typedef vector<ChromosomeMarkers> Markers;
private:
	Samples samples;
	Markers markers;
	map<string, RawSample*> byNames;
	
	RawSampleSet* clone() {
		return new RawSampleSet(*this);
	}
public:
	RawSampleSet() {
		markers.resize(nChromosomes);
	}
	RawSampleSet(SampleSet& set);
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
	RawSample* sample(const string& sampleName) {
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



class SegmentedSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class RawSampleSet;
public:
	//typedef LinearChromosome<Segment> Chromosome;
	typedef Sample<Segment, LinearChromosome<Segment> > SegmentedSample;
	
	typedef vector<SegmentedSample*> Samples;
	typedef Samples::iterator SamplesIterator;
	typedef typename SegmentedSample::Chromosomes::iterator ChromosomesIterator;
	typedef LinearChromosome<Segment>::iterator DataIterator;
private:
	Samples samples;
	map<string, SegmentedSample*> byNames;
	
	SegmentedSampleSet* clone() {
		return new SegmentedSampleSet(*this);
	}
public:
	SegmentedSampleSet() {
	}
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
	void read(const string& fileName);
	void write(const string& fileName);
	SegmentedSample* sample(const string& sampleName) {
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
};


#endif