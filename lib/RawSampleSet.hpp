#ifndef cna_RawSampleSet_h
#define cna_RawSampleSet_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "AlleleSpecific.hpp"
#include "parse.hpp"
#include "SampleSet.hpp"


// samples can be rearranged => RawSampleSet store samples as pointers
// chromosomes are not ever rearranged => Sample stores chromosomes as objects
// data can be sorted => segment data should be stored as pointers, but raw data can be stored as values


namespace cna {

class GenericSampleSet;
template <typename V> class SegmentedSampleSet;

template <typename V = rvalue>
class RawSampleSet : public SampleSet
{
	friend class GenericSampleSet;
	friend class SegmentedSampleSet<V>;
public:
	typedef SampleSet Base;
	typedef V Value;
	typedef LinearChromosome<Value> RawChromosome;
	typedef Sample<RawChromosome> RawSample;
	typedef std::vector<RawSample*> Samples;
	typedef typename Samples::iterator SamplesIterator;
	typedef typename RawSample::Chromosomes::iterator ChromosomesIterator;
	typedef typename RawChromosome::iterator DataIterator;

protected:
	Samples samples;

private:
	std::map<std::string, RawSample*> byNames;
	
	RawSampleSet* clone() const {
		return new RawSampleSet(*this);
	}
	
	void _read(std::fstream& file);
	void _write(std::fstream& file);

	void readSampleNames(FieldScanner& fields);
	void readSampleValues(FieldScanner& fields, size_t sampleStart, const std::string& chromName);

	void writeSampleNames(std::fstream& file, const char delim);
	void writeSampleValues(std::fstream& file, size_t chr, size_t markerIndex, const char delim);
	
public:
	RawSampleSet() {}
	RawSampleSet(cna::marker::Set* markerSet) : SampleSet(markerSet) {}
	RawSampleSet(const RawSampleSet& raw)
	: SampleSet(raw.markers), samples(raw.samples)  {
		byNames.clear();
		// duplicate samples
		for (size_t i = 0; i < samples.size(); ++i) {
			samples[i] = new RawSample( *(samples[i]) );
			byNames[samples[i]->name] = samples[i];
		}
		// ref the marker
		cna::marker::manager.ref(markers);
		//markers = cna::marker::manager.create(raw.markers.platform);
	}
	RawSampleSet(const cna::SegmentedSampleSet<V>& segmented);
	~RawSampleSet() {
		clear();
	}
	cna::data::Type type() {
		return cna::data::raw;
	}
	void clear() {
		SamplesIterator i, end = samples.end();
		for (i = samples.begin(); i != end; ++i) {
			// delete object pointed to by Sample* pointer
			delete (*i);
		}
		samples.clear();
		byNames.clear();
		cna::marker::manager.unref(markers);
	}
	
	RawSample* create(const std::string& sampleName) {
		RawSample* sam = byNames[sampleName];
		if (sam == NULL) {
			// sample does not exist: create it
			sam = new RawSample(sampleName);
			samples.push_back(sam);
			// register name
			byNames[sampleName] = sam;
		} else {
			// sample already exists: not supported
			throw std::runtime_error("Duplicate sample names are not permitted in RawSampleSet.");
		}
		return sam;
	}
	
	void sort();

	// TODO replace with begin() and end()
	const Samples& getSamples() const {
		return samples;
	}
	
	size_t size() {
		return samples.size();
	}
	
	void filter(const RawSampleSet& ref) {
		if (ref.markers == NULL) {
			throw std::invalid_argument("Markers in reference set are missing.");
		}
		filter(*ref.markers);
	}
	
	void filter(const cna::marker::Set& refMarkers) {
		if (markers == NULL) {
			throw std::invalid_argument("Markers in sample set are missing.");
		}
		
		// flag markers for removal
		markers->filter(refMarkers);
		
		if (samples.size() > 0) {
		
			// Create copy of samples with only unflagged markers
			
			// Shallow copy current samples (copy pointers)
			Samples oldSamples = samples;
			
			// clear current data
			samples.clear();
			byNames.clear();
			
			// pre-allocate space
			samples.resize(oldSamples.size());
			for (size_t sampleIndex = 0; sampleIndex < oldSamples.size(); ++sampleIndex) {
				samples[sampleIndex] = new RawSample(oldSamples[sampleIndex]->name);
				byNames[samples[sampleIndex]->name] = samples[sampleIndex];
				// iterate through chromosomes
				for (size_t chromIndex = 0; chromIndex < oldSamples[0]->size(); ++chromIndex) {
					// resize chromosome to be as big as chromosome from old sample
					samples[sampleIndex]->resizeChromosome(chromIndex, (*oldSamples[sampleIndex])[chromIndex].size());
				}
			}
			
			// copy data from unflagged markers
			const size_t chromEnd = markers->size();
			std::vector<size_t> validMarkersCounts(chromEnd);
			for (size_t chromIndex = 0; chromIndex < chromEnd; ++chromIndex) {
				const size_t numMarkers = (*markers)[chromIndex].size();
				validMarkersCounts[chromIndex] = 0;
				for (size_t markerIndex = 0; markerIndex < numMarkers; ++markerIndex) {
					// only copy unflagged markers
					if (!(*markers)[chromIndex][markerIndex]->flag) {
						// copy marker values from oldSamples
						const size_t numSamples = oldSamples.size();
						for (size_t sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {
							(*samples[sampleIndex])[chromIndex][validMarkersCounts[chromIndex]] = (*oldSamples[sampleIndex])[chromIndex][markerIndex];
						}
						++validMarkersCounts[chromIndex];
					}
				}
			}
			
			// shrink new sample chromosomes to match new size
			SamplesIterator it, end = samples.end();
			for (it = samples.begin(); it != end; ++it) {
				for (size_t chromIndex = 0; chromIndex < chromEnd; ++chromIndex) {
					(**it).resizeChromosome(chromIndex, validMarkersCounts[chromIndex]);
				}
			}
			
			// clear old samples
			end = oldSamples.end();
			for (it = oldSamples.begin(); it != end; ++it) {
				delete (*it);
			}
			
		}
		
		// remove flagged markers
		markers->clean();
		
	}
};

} // namespace cna

/* Template implementation */

template <typename V>
cna::RawSampleSet<V>::RawSampleSet(const cna::SegmentedSampleSet<V>& set)
{
	//TODO
	
	if (set.markers == NULL) {
		throw std::runtime_error("Marker information is required for converting from SegmentedSampleSet to RawSampleSet.");
	}
	
	throw std::runtime_error("RawSampleSet(SegmentedSampleSet<V>&) has not been implemented.");
}

template <typename V>
void cna::RawSampleSet<V>::_read(std::fstream& file)
{
	const char delim = Base::io.delim;
	const size_t nSkippedLines = Base::io.nSkippedLines, headerLine = Base::io.headerLine;
	cna::marker::Set* markers = Base::markers;
	
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	
	bool readMarkers = (markers->empty()) ? true : false;
	
	size_t lineCount = 0;
	size_t sampleStart = samples.size()-1; 
	std::string line, markerName, chromName;
	std::string_view field;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines) {
			FieldScanner fields(line, delim);
			if (lineCount == headerLine) {
				for (size_t i = 0; i < 3 && fields.next(field); ++i) {
				}
				readSampleNames(fields);
			} else {
				position pos;
				if (!fields.next(field)) continue;
				markerName.assign(field.data(), field.size());
				if (!fields.next(field)) continue;
				chromName.assign(field.data(), field.size());
				if (!fields.next(field) || !parseNumber(field, pos)) continue;
				if (readMarkers) {
					size_t chr = cna::mapping::chromosome[chromName];
					// ignore unknown chromosome: continue to next line
					if (chr == 0) continue;
					// create marker
					cna::marker::Marker* marker = new cna::marker::Marker(markerName, chr, pos);
					markers->addToChromosome(chr-1, marker);
				}
				readSampleValues(fields, sampleStart, chromName);
			}
		} else {
			// discard line
		}
	}
}

template <typename V> inline
void cna::RawSampleSet<V>::readSampleNames(FieldScanner& fields) {
	std::string_view field;
	while (fields.next(field)) {
		create(std::string(field));
	}
}

template <typename V> inline
void cna::RawSampleSet<V>::readSampleValues(FieldScanner& fields, size_t sampleStart, const std::string& chromName) {
	size_t i = sampleStart;
	Value value;
	std::string_view field;
	while (fields.next(field)) {
		if (!parseNumber(field, value)) {
			continue;
		}
		// create point at specified chromosome
		samples[++i]->addToChromosome(chromName, value);
	}
}

template <typename V>
void cna::RawSampleSet<V>::_write(std::fstream& file)
{
	const char delim = Base::io.delim;
	cna::marker::Set* markers = Base::markers;
	
	file << "marker" << delim << "chromosome" << delim << "position";
	
	writeSampleNames(file, delim);
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers->size(); ++chr) {
		for (size_t markerIndex = 0; markerIndex < markers->at(chr).size(); ++markerIndex) {
			
			// print marker information
			file << markers->at(chr)[markerIndex]->name << delim << markers->at(chr)[markerIndex]->chromosome << delim << markers->at(chr)[markerIndex]->pos;
			
			writeSampleValues(file, chr, markerIndex, delim);
		}
	}
}

template <typename V> inline
void cna::RawSampleSet<V>::writeSampleNames(std::fstream& file, const char delim) {
	// print sample names
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << std::endl;
}

template <typename V> inline
void cna::RawSampleSet<V>::writeSampleValues(std::fstream& file, size_t chr, size_t markerIndex, const char delim) {
	// iterate through samples to print values, selected the specified chromosome and marker
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it)[chr][markerIndex];
	}
	file << std::endl;
}

template <typename V>
void cna::RawSampleSet<V>::sort()
{
	cna::marker::Set* markers = Base::markers;
	
	// Construct order vector for obtaining a sorted index of markers
	cna::marker::Set::ChromosomeMarkers chromosomeMarkers;
	std::vector<RawChromosome> samplesChromosomeCopy;
	std::vector< std::pair<position, size_t> > order;
	
	for (size_t chri = 0; chri < markers->size(); ++chri) {	
		
		// Construct order vector for obtaining a sorted index of markers
		// Additionally, replicate the markers on the curent chromosome;
		//               replicate the chromosome for all samples
		
		cna::marker::Set::ChromosomeMarkers& currentMarkers = markers->at(chri);
		size_t numMarkers = currentMarkers.size();
		
		chromosomeMarkers.clear();
		samplesChromosomeCopy.clear();
		order.clear();
		
		// copy original data
		markers->copyChromosome(chri, chromosomeMarkers);
		samplesChromosomeCopy.reserve(samples.size());
		for (size_t s = 0; s < samples.size(); ++s) {
			samplesChromosomeCopy.push_back((*samples[s])[chri]);
		}
		// create order vector
		order.reserve(numMarkers);
		for (size_t j = 0; j < numMarkers; ++j) {
			order.push_back(std::make_pair(currentMarkers[j]->pos, j));
		}
		
		// Sort on the order vector instead of the original vector<Markers*>,
		//   in order to obtain the sorted index
		std::sort(order.begin(), order.end(), &cna::compare::pair<position, size_t>);
		
		// Now, order is sorted by chromosome position, and order[i].second
		//   contains each sorted index
		
		// Sort the markers on current chromosome, and each sample
		for (size_t j = 0; j < numMarkers; ++j) {
			size_t index = order[j].second;
			// Set the marker to the corresponding sorted marker
			currentMarkers[j] = chromosomeMarkers[index];
			
			// Iterate through samples, set the corresponding marker in the current chromosome
			for (size_t s = 0; s < samples.size(); ++s) {
				(*samples[s])[chri][j] = samplesChromosomeCopy[s][index];
			}
		}
		
	}
}


/* Template Specialization */

// Use the same specialization for alleles_cn and alleles_rcn

#define SPECIALIZATION_TYPE alleles_cn
#include "RawSampleSet_special.hpp"
#undef SPECIALIZATION_TYPE

#define SPECIALIZATION_TYPE alleles_rcn
#include "RawSampleSet_special.hpp"
#undef SPECIALIZATION_TYPE


#endif
