#ifndef genomic_Sample_h
#define genomic_Sample_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "typedefs.h"
#include "global.hpp"
#include "Chromosome.hpp"
#include "CopyNumberValue.hpp"
#include "Segment.hpp"

template <typename Chromosome>
class Sample
{
public:
	typedef typename Chromosome::DataType T;
	typedef std::vector<Chromosome> Chromosomes;
	std::string name;
	std::string platform;
private:
	// vector of segments vectors for each Chromosome
	Chromosomes items;
public:
	
	Sample(const std::string& sampleName) : name(sampleName) {
		// always allocate for all chromosomes
		items.reserve(nChromosomes);
		for (chromid i = 0; i < nChromosomes; ++i) {
			items.push_back(Chromosome(i));
		}
	}
	
	//Sample(const Sample& sample)
	//: name(sample.name), platform(sample.platform), items(sample.items) {}
	
	~Sample() {
		clear();
	}
	void clear() {
		items.clear();
	}
	const size_t size() const {
		items.size();
	}
	Chromosome& at(chromid chromIndex) {
		return this->operator[](chromIndex);
	}
	Chromosome& at(const std::string& chromName) {
		return this->operator[](chromName);
	}
	Chromosome& operator[](chromid chromIndex) {
		if (chromIndex < items.size()) {
			return items[chromIndex];
		} else {
			throw logic_error("Chromosome index is out of bound.");
		}
	}
	Chromosome& operator[](const std::string& chromName) {
		chromid k = mapping::chromosome[chromName];
		if (k != 0) return items[k-1];
		throw logic_error("Chromosome name is not found.");
	}
	typename Chromosomes::iterator begin() {
		return items.begin();
	}
	typename Chromosomes::iterator end() {
		return items.end();
	}
	typename Chromosomes::const_iterator begin() const {
		return items.begin();
	}
	typename Chromosomes::const_iterator end() const {
		return items.end();
	}
	Chromosome* chromosome(chromid chromIndex) {
		if (chromIndex < items.size()) {
			return &(items[chromIndex]);
		}
		return NULL;
	}
	Chromosome* chromosome(const std::string& chromName) {
		chromid k = mapping::chromosome[chromName];
		if (k != 0) return &(items[k-1]);
		return NULL;
	}
	void addToChromosome(chromid chromIndex, const T& item) {
		if (chromIndex < items.size()) {
			items[chromIndex].push_back(item);
		}
	}
	void addToChromosome(const std::string& chromName, const T& item) {
		chromid k = mapping::chromosome[chromName];
		if (k != 0) {
			items[k-1].push_back(item);
		}
	}
	void resizeChromosome(chromid chromIndex, size_t n) {
		if (chromIndex < items.size()) {
			items[chromIndex].resize(n);
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


#endif
