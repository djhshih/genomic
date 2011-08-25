#ifndef genomic_GenomicGraph_h
#define genomic_GenomicGraph_h

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "Graph.hpp"


// handles RawSampleSet and SegmentedSampleSet

class GenomicGraph : Graph
{
public:
	
	template <typename V>
	void add(const RawSampleSet<V>& set) {
		// data sanity check
		
		plots.push_back( boost::bind(&GenomicGraph::plot, this, set) );
	}
	
	template <typename V>
	void plot(const RawSampleSet<V>& set) {
		
	}
	
	
	template <typename V>
	void add(const SegmentedSampleSet<V>& set) {
		// data sanity check
	
		plots.push_back( boost::bind(&GenomicGraph::plot, this, set) );
	}
	
	template <typename V>
	void plot(const SegmentedSampleSet<V>& set) {
		
	}
	
private:	
	
	
	
}

#endif
