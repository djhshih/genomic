#ifndef genomic_Properties_h
#define genomic_Properties_h

struct IOProperties
{
	char delim;
	size_t nSkippedLines;
	size_t headerLine;
	bool mergeSamples;
	
	IOProperties(char _delim='\t',
	             size_t _headerLine=1,
	             size_t _nSkippedLines=0,
	             bool _mergeSamples=false)
	: delim(_delim),
	  headerLine(_headerLine),
	  nSkippedLines(_nSkippedLines),
	  mergeSamples(_mergeSamples)
	{}
	
};

struct CNACriteria
{
	float reference;
	float deviation;
	
	CNACriteria() : reference(-1), deviation(-1) {}
	
	CNACriteria(float referenceState, float stateDeviation)
	: reference(referenceState), deviation(stateDeviation)
	{}
};

#endif
