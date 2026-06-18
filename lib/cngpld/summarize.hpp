#ifndef cngpld_summarize_h
#define cngpld_summarize_h

#include <string>
#include <vector>

#include "typedefs.h"
#include "SegmentedSampleSet.hpp"

namespace cngpld {

struct CNSummaryPoint {
	position pos;
	double value;
};

typedef std::vector<CNSummaryPoint> CNSummary;

double summarize_cn_at_position(
	const cna::SegmentedSampleSet<rvalue>& seg,
	const std::string& sample,
	size_t chrom_index,
	position pos,
	int direction,
	double cutoff);

CNSummary summarize_cn(
	const cna::SegmentedSampleSet<rvalue>& seg,
	const std::string& sample,
	size_t chrom_index,
	int direction,
	double cutoff,
	const std::vector<position>* positions = NULL);

} // namespace cngpld

#endif
