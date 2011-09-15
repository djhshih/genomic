#ifndef genomic_SegmentedSampleSet_h
#define genomic_SegmentedSampleSet_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "AlleleSpecific.hpp"
#include "SampleSet.hpp"

template <typename V> class RawSampleSet;

template <typename V = rvalue>
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
	std::map<string, SegmentedSample*> byNames;
	
	SegmentedSampleSet* clone() const {
		return new SegmentedSampleSet(*this);
	}
	
	void _read(fstream& file);
	void _write(fstream& file);
	
	void readSegment(istringstream& stream, Segment<V>& seg) {
		stream >> seg.start >> seg.end;
		if (positionsOnly) {
			seg.count = seg.value = 0;
		} else {
			stream >> seg.count >> seg.value;
		}
	}
	
	size_t _find(SegmentedChromosome& array, position x);
	
public:
	SegmentedSampleSet() {
		_setIO();
	}
	SegmentedSampleSet(marker::Set* markerSet) : SampleSet(markerSet) {
		_setIO();
	}
	SegmentedSampleSet(const SegmentedSampleSet& segmented)
	: SampleSet(segmented.markers), samples(segmented.samples) {
		_setIO();
		byNames.clear();
		// duplicate samples
		for (size_t i = 0; i < samples.size(); ++i) {
			samples[i] = new SegmentedSample( *(samples[i]) );
			byNames[samples[i]->name] = samples[i];
		}
		// ref the marker
		marker::manager.ref(markers);
		//markers = marker::manager.create(raw.markers.platform);
	}
	SegmentedSampleSet(const RawSampleSet<V>& raw);
	~SegmentedSampleSet() {
		clear();
	}
	data::Type type() {
		return data::segmented;
	}
	void clear() {
		typename Samples::iterator it, end = samples.end();
		for (it = samples.begin(); it != end; ++it) {
			// delete object pointed to by Sample* pointer
			delete (*it);
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
	
	size_t size() {
		return samples.size();
	}
	
	size_t find(const string& sampleName, size_t chromIndex, position start) {
		return _find(*(byNames[sampleName]->chromosome(chromIndex)), start);
	}
	size_t find(SegmentedSample* sample, size_t chromIndex, position start) {
		return _find(*(sample->chromosome(chromIndex)), start);
	}
	
	void markAberrant(float refValue, float diff);
	
	void filter(SegmentedSampleSet& ref, float diceThreshold, bool merge, bool aberrantOnly);
	
protected:
	bool mergeSamples;
	bool positionsOnly;
	
private:
	void _setIO() {
		mergeSamples = false;
		positionsOnly = false;
	}
	
};


/* Template implementation */

template <typename V>
SegmentedSampleSet<V>::SegmentedSampleSet(const RawSampleSet<V>& raw)
{
	clear();
	// use iterators to avoid assuming RawSampleSet stores data in vectors
	// however, need to assume that markers are stored in vectors, for looking up marker information
	
	typedef typename RawSampleSet<V>::Samples::const_iterator RawSamplesIterator;
	typedef typename RawSampleSet<V>::RawSample::Chromosomes::const_iterator RawChromosomesIterator;
	typedef typename RawSampleSet<V>::RawChromosome::const_iterator RawDataIterator;
	
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
							raw.markers->at(chr)[startMarkerIndex]->pos,
							raw.markers->at(chr)[markerIndex-1]->pos,
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
					raw.markers->at(chr)[startMarkerIndex]->pos,
					raw.markers->at(chr)[markerIndex-1]->pos,
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
	//const char delim = Base::delim;
	//const size_t nSkippedLines = Base::nSkippedLines, headerLine = Base::headerLine;
	
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	
	size_t lineCount = 0;
	string line, sampleName, chromName;
	Segment<Value>* seg;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > io.nSkippedLines && lineCount != io.headerLine) {
			istringstream stream(line);
			stream >> sampleName >> chromName;
			// ignore unknown chromosome: continue to next line
			if (mapping::chromosome[chromName] == 0) continue;
			// create segment at specified chromosome
			Segment<V> seg;
			readSegment(stream, seg);
			if (mergeSamples) sampleName = "ALL";
			create(sampleName)->addToChromosome(chromName, seg);
			//trace("%s %s %d %d %d %f\n", sampleName.c_str(), chromName.c_str(), seg.start, seg.end, seg.count, seg.value);
		} else {
			// discard line
		}
	}
}

template <typename V>
void SegmentedSampleSet<V>::_write(fstream& file)
{
	const char delim = Base::io.delim;
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "state" << endl;
	
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->count << delim << segIt->value << endl;
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
void SegmentedSampleSet<V>::markAberrant(float refValue, float diff) {
	size_t count =  0;
	// iterate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				if (segIt->value <= refValue - diff || segIt->value >= refValue + diff) {
					segIt->aberrant = true;
					++count;
				} else {
					segIt->aberrant = false;
				}
			}
		}
	}
	trace("Number of segments marked aberrant: %d\n", count);
}

template <typename V>
void SegmentedSampleSet<V>::filter(SegmentedSampleSet& ref, float diceThreshold, bool merge=false, bool aberrantOnly=false)
{
	size_t filteredCount = 0;
	Samples oldSamples;
	// iterate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				
				if (!aberrantOnly || segIt->aberrant) {
					// Compare against segments in all samples in reference set
					SamplesIterator refIt, refEnd = ref.samples.end();
					bool filterSegment = false;
					for (refIt = ref.samples.begin(); refIt != refEnd; ++refIt) {
						
						// determine lower and upper bounds
						SegmentedChromosome& refChrom = (**refIt)[chri];
						position_diff lower = 2*(diceThreshold-1)/diceThreshold*segIt->end + (2-diceThreshold)/diceThreshold*(segIt->start - 1) + 1;
						if (lower < 0) lower = 0;
						position_diff upper = 2*(1-diceThreshold)/(2-diceThreshold)*segIt->end + diceThreshold/(2-diceThreshold)*(segIt->start - 1) + 1;
						//cout << segIt->start << " " << segIt->end << " " << lower << " " << upper << endl;
						
						// find marker indices corresponding to lower and upper bounds
						size_t lowerIndex = ref.find(*refIt, chri, lower);
						size_t upperIndex = ref.find(*refIt, chri, upper) + 1;
						if (upperIndex >= refChrom.size()) upperIndex = refChrom.size()-1;
						
						//size_t lowerIndex = 0, upperIndex = refChrom.size()-1;
						//cout << "Index: " << lowerIndex << ", " << upperIndex << endl;
						for (size_t i = lowerIndex; i <= upperIndex; ++i) {
							// calculate Dice coefficient
							position_diff intersection = min(refChrom[i].end, segIt->end) - max(refChrom[i].start, segIt->start) + 1;
							//cout << segIt->start << " " << refChrom->at(i).start << " " << intersection << endl;
							if (intersection > 0) {
								float dice = 2 * float(intersection) / (refChrom[i].length() + segIt->length());
								if (dice > diceThreshold) {
									//cout << "Filter: " << segIt->start << " " << refChrom->at(i).start << " " << dice << endl;
									// Mark segment for deletion
									segIt->flag = true;
									filterSegment = true;
									trace("Filter chr%s:%d-%d in %s\n", mapping::chromosome[chri+1].c_str(), segIt->start, segIt->end, (*it)->name.c_str());
									++filteredCount;
									break;
								}
							}
						}
						if (filterSegment) break;
					}
				}  // if (aberrantOnly && segIt->aberrant)
				
			}
			++chri;
		}
		
		// copy current samples
		oldSamples.push_back(*it);
	}
	
	// Create new sample set with marked segments removed or merged
	samples.clear();
	byNames.clear();
	end = oldSamples.end();
	for (it = oldSamples.begin(); it != end; ++it) {
		// create new sample
		SegmentedSample* sample = create((*it)->name);
		
		ChromosomesIterator chrIt;
		const ChromosomesIterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			DataIterator segIt;
			const DataIterator segEnd = chrIt->end();
			string chrom = mapping::chromosome[chri+1];
			
			Segment<V>* prevUnmarkedSegment = NULL, *nextUnmarkedSegment;
			// prevUnmarkedSegment will point to previous unmarked segment in the samples
			//  so that the unmarked segment is modified after being added to samples
			// nextUnmarkedSegment will point to the next unmarked segment in the oldSamples
			//  so that the unmarked segment is modified before being added to samples
			
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				// only copy unflagged segments
				if (!segIt->flag) {
					Segment<V> seg(segIt->start, segIt->end, segIt->count, segIt->value);
					prevUnmarkedSegment = sample->addToChromosome(chri, seg);
					//prevUnmarkedSegment = &(*segIt);
				} else if (merge) {
					// segment is flagged for deletion
					// option to merge is enabled: merge segment with upstream or downstream
					//  segment, whichever is bigger
					
					// find next unmarked segment
					nextUnmarkedSegment = NULL;
					DataIterator tmp = segIt;
					do {
						if (++tmp != segEnd) {
							if (!tmp->flag) {
								nextUnmarkedSegment = &(*tmp);
								break;
							}
						} else {
							break;
						}
					} while (nextUnmarkedSegment == NULL);
					
					if (prevUnmarkedSegment == NULL && nextUnmarkedSegment == NULL) {
						trace("Warning: segment chr%s:%d-%d in %s cannot be merge with another segment\n",
								chrom.c_str(), segIt->start, segIt->end, (*it)->name.c_str());
					} else {
						
						// Compare upstream and downstream unmarked segments
						// skip adding 1 to get correct size
						position_diff prevSize, nextSize;
						
						if (prevUnmarkedSegment != NULL) {
							prevSize = prevUnmarkedSegment->end - prevUnmarkedSegment->start;
						} else {
							prevSize = 0;
						}
						
						if (nextUnmarkedSegment != NULL) {
							nextSize = nextUnmarkedSegment->end - nextUnmarkedSegment->start;
						} else {
							nextSize = 0;
						}
						
						if (prevSize >= nextSize) {
							
							trace("Merge chr%s:%d-%d to upstream chr%s:%d-%d in %s\n",
								chrom.c_str(), segIt->start, segIt->end,
								chrom.c_str(), prevUnmarkedSegment->start, prevUnmarkedSegment->end,
								(*it)->name.c_str() );
							
							// extend upstream segment
							prevUnmarkedSegment->end = segIt->end;
							
						} else {
							// extend downstream segment
							if (nextUnmarkedSegment->start > segIt->start) {
								// check guards against multiple assignments in cases where a series
								//  of segments are marked for deletion
								
								trace("Merge chr%s:%d-%d to downstream chr%s:%d-%d in %s\n",
									chrom.c_str(), segIt->start, segIt->end,
									chrom.c_str(), nextUnmarkedSegment->start, nextUnmarkedSegment->end,
									(*it)->name.c_str() );
								
								nextUnmarkedSegment->start = segIt->start;
							}
						}
						
					}
				} // if (!segIt->flag)
			}
			++chri;
		}
		// delete old sample
		delete (*it);
	}
	oldSamples.clear();
	
	trace("Number of segments filtered: %d\n", filteredCount);
}


/* Template Specialization */

// Use the same specialization for alleles_cn and alleles_rcn

#define SPECIALIZATION_TYPE alleles_cn
#include "SegmentedSampleSet_special.hpp"
#undef SPECIALIZATION_TYPE

#define SPECIALIZATION_TYPE alleles_rcn
#include "SegmentedSampleSet_special.hpp"
#undef SPECIALIZATION_TYPE

#endif
