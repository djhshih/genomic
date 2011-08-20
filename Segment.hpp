#ifndef genomic_Segment_h
#define genomic_Segment_h

#include "typedefs.h"

template <typename V>
class Segment
{
public:
	position start;
	position end;
	position count;
	V value;
	bool flag;
	position length() {
		return end - start + 1;
	}
	Segment() : flag(false) {}
	Segment(position startPos, position endPos, unsigned long numElements, V segValue) : flag(false), start(startPos), end(endPos), count(numElements), value(segValue) {}
	static bool compare(const Segment& a, const Segment& b) {
		return a.start < b.start; 
	}
};

#endif
