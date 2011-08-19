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

#include "global.hpp"
#include "Chromosome.hpp"
#include "CopyNumberValue.hpp"
#include "Segment.hpp"


//TODO use boost::spirit for parsing
//TODO allow use of input delimiter other than whitespace in all functions
//TODO improve sorting (enable skipping): used up 50GB of RAM on 3 picnic samples!


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
	Sample(const std::string& sampleName) : name(sampleName), items(nChromosomes) {
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
	Chromosome* operator[](const std::string& chromName) {
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
	Chromosome* chromosome(const std::string& chromName) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) return &(items[k-1]);
		return NULL;
	}
	void addToChromosome(size_t chromIndex, const T& item) {
		if (chromIndex < items.size()) {
			items[chromIndex].push_back(item);
		}
	}
	void addToChromosome(const std::string& chromName, const T& item) {
		size_t k = mapping::chromosome[chromName];
		if (k != 0) {
			items[k-1].push_back(item);
		}
	}
	void resizeChromosome(size_t chromIndex, size_t n) {
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
