#ifndef genomic_RawSampleSet_h
#define genomic_RawSampleSet_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "CopyNumberValue.hpp"
#include "SampleSet.hpp"


// samples can be rearranged => RawSampleSet store samples as pointers
// chromosomes are not ever rearranged => Sample stores chromosomes as objects
// data can be sorted => store segment data as pointers, but raw data can be stored as values

extern marker::Manager marker::manager;

template <typename T> class LinearChromosome;
template <typename Chromosome> class Sample;
class GenericSampleSet;
template <typename V> class SegmentedSampleSet;

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
	std::map<string, RawSample*> byNames;
	
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
	size_t size() {
		return samples.size();
	}
	
	void filter(const RawSampleSet& ref) {
		if (ref.markers == NULL) {
			throw logic_error("Markers in reference set are missing");
		}
		filter(*ref.markers);
	}
	
	void filter(const marker::Set& refMarkers) {
		if (markers == NULL) {
			throw logic_error("Markers in sample set are missing.");
		}
		
		// flag markers for removal
		markers->filter(refMarkers);
		
		if (samples.size() > 0) {
		
			// Create copy of samples with only unflagged markers
			
			// Shallow copy current samples (copy pointers)
			Samples oldSamples = samples;
			
			samples.clear();
			byNames.clear();
			
			samples.resize(oldSamples.size());
			for (size_t sampleIndex = 0; sampleIndex < oldSamples.size(); ++sampleIndex) {
				// iterate through chromosomes
				for (size_t chromIndex = 0; chromIndex < oldSamples[0]->size(); ++chromIndex) {
					// resize chromosome to be as big as chromosome from old sample
					samples[sampleIndex]->resizeChromosome(chromIndex, (*oldSamples[sampleIndex])[chromIndex].size());
				}
			}
			
			const size_t chromEnd = markers->size();
			vector<size_t> validMarkersCounts(chromEnd);
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
			
			// resize new sample chromosomes
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



/* Template implementation */

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
					marker::Marker* marker = new marker::Marker(markerName, chr, pos);
					markers->addToChromosome(chr-1, marker);
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
			file << markers->at(chr)[markerIndex]->name << delim << markers->at(chr)[markerIndex]->chromosome << delim << markers->at(chr)[markerIndex]->pos;
			
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
		file << delim << (**it)[chr][markerIndex];
	}
	file << endl;
}

// TODO make more efficient!
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
			order.push_back( make_pair(markers->at(chri)[j]->pos, j) );
			
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
				(*samples[s])[chri][j] = samplesChromosomeCopy[s][index];
			}
		}
		
	}
}


/* Template Specialization */

// Use the same specialization for AlleleSpecificCopyNumberValue and AlleleSpecificIntegerCopyNumberValue

#define SPECIALIZATION_TYPE AlleleSpecificCopyNumberValue
#include "RawSampleSet.special"
#undef SPECIALIZATION_TYPE

#define SPECIALIZATION_TYPE AlleleSpecificIntegerCopyNumberValue
#include "RawSampleSet.special"
#undef SPECIALIZATION_TYPE


#endif
