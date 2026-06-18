#include "cngpld/summarize.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

const cna::SegmentedSampleSet<rvalue>::Segments& get_segments(
	const cna::SegmentedSampleSet<rvalue>& segset,
	const std::string& sample_name,
	size_t chrom_index)
{
	const cna::SegmentedSampleSet<rvalue>::SegmentedSample* sam = segset.sample(sample_name);
	if (sam == NULL) {
		throw std::invalid_argument("Unknown sample: '" + sample_name + "'.");
	}
	if (chrom_index >= sam->size()) {
		throw std::invalid_argument("Chromosome index is out of bounds.");
	}
	return sam->at(static_cast<chromid>(chrom_index));
}

std::vector<position> default_positions(const cna::SegmentedSampleSet<rvalue>::Segments& segments)
{
	std::vector<position> positions;
	positions.reserve(segments.size() * 2);
	for (cna::SegmentedSampleSet<rvalue>::Segments::const_iterator it = segments.begin(); it != segments.end(); ++it) {
		positions.push_back(it->start);
		positions.push_back(it->end);
	}
	std::sort(positions.begin(), positions.end());
	positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
	return positions;
}

}

namespace cngpld {

double summarize_cn_at_position(
	const cna::SegmentedSampleSet<rvalue>& seg,
	const std::string& sample,
	size_t chrom_index,
	position pos,
	int direction,
	double cutoff)
{
	if (direction != 1 && direction != -1) {
		throw std::invalid_argument("direction must be 1 or -1.");
	}

	const cna::SegmentedSampleSet<rvalue>::Segments& segments = get_segments(seg, sample, chrom_index);
	size_t overlap_count = 0;
	size_t altered_count = 0;
	double altered_sum = 0.0;

	for (cna::SegmentedSampleSet<rvalue>::Segments::const_iterator it = segments.begin(); it != segments.end(); ++it) {
		if (it->start > it->end) {
			throw std::invalid_argument("Segment start is greater than end.");
		}
		if (it->start <= pos && pos <= it->end) {
			++overlap_count;
			const double adj = static_cast<double>(direction) * static_cast<double>(it->value);
			if (adj > cutoff) {
				altered_sum += std::exp(adj);
				++altered_count;
			}
		}
	}

	if (altered_count == 0) {
		return 0.0;
	}
	return altered_sum / static_cast<double>(overlap_count);
}

CNSummary summarize_cn(
	const cna::SegmentedSampleSet<rvalue>& seg,
	const std::string& sample,
	size_t chrom_index,
	int direction,
	double cutoff,
	const std::vector<position>* positions)
{
	const cna::SegmentedSampleSet<rvalue>::Segments& segments = get_segments(seg, sample, chrom_index);
	std::vector<position> computed_positions;
	if (positions == NULL) {
		computed_positions = default_positions(segments);
		positions = &computed_positions;
	}

	CNSummary out;
	out.reserve(positions->size());
	for (std::vector<position>::const_iterator it = positions->begin(); it != positions->end(); ++it) {
		CNSummaryPoint pt;
		pt.pos = *it;
		pt.value = summarize_cn_at_position(seg, sample, chrom_index, *it, direction, cutoff);
		out.push_back(pt);
	}
	return out;
}

} // namespace cngpld
