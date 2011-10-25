#ifndef genomic_Segment_h
#define genomic_Segment_h

#include "typedefs.h"

template <typename V>
class Segment
{
public:
	// chromosome number (1-index)
	chromid chromosome;
	position start;
	position end;
	position count;
	V value;
	bool flag;
	bool aberrant;
	bool valid;
	
	position length() {
		return end - start + 1;
	}
	
	Segment() : flag(false), aberrant(false), valid(true) {}
	
	Segment(chromid chromosomeNumber)
	: flag(false), aberrant(false), valid(true),
	chromosome(chromosomeNumber) {}
	
	Segment(chromid chromosomeNumber, position startPos, position endPos, unsigned long numElements, V segValue)
	: flag(false), aberrant(false), valid(true),
	  chromosome(chromosomeNumber),
	  start(startPos), end(endPos), count(numElements), value(segValue) {}
	  
	static bool compare(const Segment& a, const Segment& b) {
		return a.start < b.start; 
	}
};

#endif
