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

 
template <typename V>
class filter_operator
{
public:
	filter_operator() {}
	bool operator()(const Segment<V>& seg) {
		return false;
	}
};

template <typename V>
class spurious_segment_filter : filter_operator<V>
{
	position count;
public:
	spurious_segment_filter(position countThreshold) : count(countThreshold) {}
	bool operator()(const Segment<V>& seg) {
		return seg.count < count;
	}
};

template <typename V>
class small_segment_filter : filter_operator<V>
{
	position length;
public:
	small_segment_filter(size_t lengthThreshold) : length(lengthThreshold) {}
	bool operator()(const Segment<V>& seg) {
		return (seg.end - seg.start + 1) < length;
	}
};

template <typename V, typename T=float>
class balanced_segment_filter : filter_operator<V>
{
	// note: T cannot be an unsigned type; operator+ and operator- must be defined
	// TODO: use one template type instead of two
	T reference;
	T deviation;
public:
	balanced_segment_filter(T referenceState, T stateDeviation)
	: reference(referenceState), deviation(stateDeviation) {}
	bool operator()(const Segment<V>& seg) {
		return !(seg.value <= reference - deviation && seg.value >= reference + deviation);
	}
};

template <typename V = rvalue>
class SegmentedSampleSet : public SampleSet
{
	
	friend class GenericSampleSet;
	friend class RawSampleSet<V>;
	
public:
	
	typedef SampleSet Base;
	typedef V Value;
	typedef LinearChromosome< Segment<Value> > Segments;
	typedef Sample<Segments> SegmentedSample;
	
	typedef vector<SegmentedSample*> Samples;
	typedef typename SegmentedSample::Chromosomes Chromosomes;
	
// 	typedef typename Samples::iterator SamplesIterator;
// 	typedef typename Chromosomes::iterator ChromosomesIterator;
// 	typedef typename Segments::iterator DataIterator;
	
	
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
	
	void markAberrant();
	
	void removeFlagged(bool merged);
	
	size_t _find(Segments& array, position x);
	
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
	
	void set(const CNACriteria& criteria) {
		cna = criteria;
	}
	
	void reset();
	
	template <typename filter_operator_type>
	void filter(filter_operator_type f, bool merge);
	
	void filter(SegmentedSampleSet& ref, float diceThreshold, bool merge, bool aberrantOnly, bool optimize);
	
	
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
	
	typename Samples::iterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		typename Chromosomes::iterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			typename Segments::iterator segIt, segEnd = chrIt->end();
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
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::iterator chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			std::sort(chrIt->begin(), chrIt->end(), &Segment<Value>::compare);
		}
	}
}

// find segments that start at specified position
template <typename V>
size_t SegmentedSampleSet<V>::_find(Segments& array, position x) {
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
void SegmentedSampleSet<V>::markAberrant() {
	size_t count =  0;
	// iterate through samples
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::const_iterator chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			typename Segments::iterator segIt;
			typename Segments::const_iterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				if (segIt->value <= cna.reference - cna.deviation || segIt->value >= cna.reference + cna.deviation) {
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
void SegmentedSampleSet<V>::reset() {
	// iterate through samples
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::const_iterator chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			typename Segments::iterator segIt;
			typename Segments::const_iterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				segIt->flag = false;
				segIt->aberrant = false;
				segIt->valid = true;
			}
		}
	}
}

template <typename V>
template <typename filter_operator_type>
void SegmentedSampleSet<V>::filter(filter_operator_type f, bool merge=false)
{
	size_t filteredCount = 0;
	// iterate through samples
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::const_iterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			typename Segments::iterator segIt;
			typename Segments::const_iterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				// use filter functor to flag segment
				segIt->flag = f(*segIt);
				if (segIt->flag) ++filteredCount;
			}
		}
	}
	
	removeFlagged(merge);
	trace("Number of segments filtered: %d\n", filteredCount);
}

// assume flags of all segments are set to false
template <typename V>
void SegmentedSampleSet<V>::filter(SegmentedSampleSet& ref, float diceThreshold, bool merge=false, bool aberrantOnly=false, bool optimize=true)
{
	if (aberrantOnly) markAberrant();
	
	size_t filteredCount = 0;
	// iterate through samples
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::const_iterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			const char* chrom = mapping::chromosome[chri+1].c_str();
			
			// iterate through segments on a chromosome
			typename Segments::iterator segIt;
			typename Segments::const_iterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				
				if (!aberrantOnly || segIt->aberrant) {
					// Compare against segments in all samples in reference set
					typename Samples::iterator refIt;
					typename Samples::const_iterator refEnd = ref.samples.end();
					bool filterSegment = false;
					for (refIt = ref.samples.begin(); refIt != refEnd; ++refIt) {
						
						// determine lower and upper bounds
						Segments& refChrom = (**refIt)[chri];
						if (refChrom.size() > 0) {
							// chromosome may be empty
							
							size_t lowerIndex, upperIndex;
							
							if (optimize) {
								position_diff lower = 2*(diceThreshold-1)/diceThreshold*segIt->end + (2-diceThreshold)/diceThreshold*(segIt->start - 1) + 1;
								if (lower < 0) lower = 0;
								position_diff upper = 2*(1-diceThreshold)/(2-diceThreshold)*segIt->end + diceThreshold/(2-diceThreshold)*(segIt->start - 1) + 1;
								//cout << segIt->start << " " << segIt->end << " " << lower << " " << upper << endl;
								
								// find marker indices corresponding to lower and upper bounds
								lowerIndex = ref.find(*refIt, chri, lower);
								upperIndex = ref.find(*refIt, chri, upper) + 1;
								if (upperIndex >= refChrom.size()) upperIndex = refChrom.size()-1;
							} else {
								lowerIndex = 0;
								upperIndex = refChrom.size()-1;
							}
							
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
										trace("Filter chr%s:%d-%d in %s: %.2f overlap with chr%s:%d-%d in reference\n",
													chrom, segIt->start, segIt->end, (*it)->name.c_str(),
													dice,
													chrom, refChrom[i].start, refChrom[i].end);
										++filteredCount;
										break;
									}
								}
							}
						}
						if (filterSegment) break;
					}
				}  // if (aberrantOnly && segIt->aberrant)
				
			}
			++chri;
		}
	}
	
	removeFlagged(merge);
	
	trace("Number of segments filtered: %d\n", filteredCount);
}

template <typename V>
void SegmentedSampleSet<V>::removeFlagged(bool merge)
{
	Samples oldSamples;
	
	// copy current samples (pointers)
	typename Samples::iterator it;
	typename Samples::const_iterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		oldSamples.push_back(*it);
	}
	
	// clear vectors, but do not deallocate, since oldSamples hold pointers to same objects as samples
	samples.clear();
	byNames.clear();
	
	// Create new sample set with marked segments removed or merged
	end = oldSamples.end();
	for (it = oldSamples.begin(); it != end; ++it) {
		// create new sample
		SegmentedSample* sample = create((*it)->name);
		
		typename Chromosomes::iterator chrIt;
		typename Chromosomes::const_iterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			typename Segments::iterator segIt;
			typename Segments::const_iterator segEnd = chrIt->end();
			const char* chrom = mapping::chromosome[chri+1].c_str();
			
			Segment<V>* prevUnmarkedSegment = NULL, *nextUnmarkedSegment;
			// prevUnmarkedSegment will point to previous unmarked segment in the samples
			//  so that the unmarked segment is modified after being copied to samples
			// nextUnmarkedSegment will point to the next unmarked segment in the oldSamples
			//  so that the unmarked segment is modified before being copied to samples
			
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				
				if (!segIt->flag) {
					
					// only create new copy of unflagged segments
					Segment<V> seg(segIt->start, segIt->end, segIt->count, segIt->value);
					prevUnmarkedSegment = sample->addToChromosome(chri, seg);
					
				} else if ( segIt->valid && merge ) {
					// segment is flagged for deletion
					// option to merge is enabled: merge segment with upstream or downstream
					//  segment, whichever is bigger
					
					// only valid segments are candidates for merging, as guard against merging of consecutive segments,
					//   which results in undesirable behaviour
					// since only flagged segments are ever marked as valid
					// copying only unflagged segments ensure that all segments are valid
					
					// find next unmarked segment
					nextUnmarkedSegment = NULL;
					typename Segments::iterator tmp = segIt;
					do {
						if (++tmp != segEnd) {
							if (tmp->flag) {
								// mark traversed flagged segments as invalid
								//  s.t. they are not subsequent candidate for merging
								tmp->valid = false;
							} else {
								// segment is unmarked
								nextUnmarkedSegment = &(*tmp);
								break;
							}
						} else {
							break;
						}
					} while (nextUnmarkedSegment == NULL);
					
					if (prevUnmarkedSegment == NULL && nextUnmarkedSegment == NULL) {
						
						trace("Warning: segment chr%s:%d-%d in %s cannot be merge with another segment\n",
								chrom, segIt->start, segIt->end, (*it)->name.c_str());
						
					} else {
						
						if ( prevUnmarkedSegment != NULL && nextUnmarkedSegment != NULL &&
							   absdiff(prevUnmarkedSegment->value, nextUnmarkedSegment->value) <= cna.deviation ) {
							// previous and next segments are essentially the same
							
							trace("Merge segments from chr%s:%d-%d to chr%s:%d-%d in %s\n",
								chrom, prevUnmarkedSegment->start, prevUnmarkedSegment->end,
								chrom, nextUnmarkedSegment->start, nextUnmarkedSegment->end,
								(*it)->name.c_str() );
							
							// merge next unmarked segment to previous unmarked segment
							prevUnmarkedSegment->end = nextUnmarkedSegment->end;
							prevUnmarkedSegment->count += nextUnmarkedSegment->count;
							float totalCount = prevUnmarkedSegment->count + nextUnmarkedSegment->count;
							
							// update value with weighted average
							if (neq(prevUnmarkedSegment->value, nextUnmarkedSegment->value)) {
								prevUnmarkedSegment->value = 
									prevUnmarkedSegment->value * (prevUnmarkedSegment->count/totalCount) + 
									nextUnmarkedSegment->value * (nextUnmarkedSegment->count/totalCount);
							}
							
							// mark the segment for removal, since it has been merged
							// also mark it as invalid, s.t. it is not a subsequent candidate for merging
							nextUnmarkedSegment->flag = true;
							nextUnmarkedSegment->valid = false;
							
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
									chrom, segIt->start, segIt->end,
									chrom, prevUnmarkedSegment->start, prevUnmarkedSegment->end,
									(*it)->name.c_str() );
								
								// extend upstream segment
								prevUnmarkedSegment->end = segIt->end;
								// do not increase segment count, because that'd be lying
								
							} else {
								if (nextUnmarkedSegment->start > segIt->start) {
									// check guards against multiple assignments in cases where a series
									//  of segments are marked for deletion
									
									trace("Merge chr%s:%d-%d to downstream chr%s:%d-%d in %s\n",
										chrom, segIt->start, segIt->end,
										chrom, nextUnmarkedSegment->start, nextUnmarkedSegment->end,
										(*it)->name.c_str() );
									
									// extend downstream segment
									nextUnmarkedSegment->start = segIt->start;
									// do not increase segment count, because that'd be lying
								}
							}
							
						}  // if (prevUnmarkedSegment != NULL & nextUnmarkedSegment != NULL)
						
					}
				} // if (!segIt->flag)
			}
			++chri;
		}
		// delete old sample
		delete (*it);
	}
	oldSamples.clear();
	
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
