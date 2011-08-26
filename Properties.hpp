#ifndef genomic_Properties_h
#define genomic_Properties_h

struct PropertiesIO {
	
	char delim;
	size_t nSkippedLines;
	size_t headerLine;
	bool mergeSamples;
	
	PropertiesIO(char _delim='\t',
	             size_t _headerLine=1,
	             size_t _nSkippedLines=0,
	             bool _mergeSamples=false)
	: delim(_delim),
	  headerLine(_headerLine),
	  nSkippedLines(_nSkippedLines),
	  mergeSamples(_mergeSamples)
	{}
	
};

#endif
