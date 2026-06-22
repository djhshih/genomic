#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Unit Tests for cngpld"
#include <boost/test/unit_test.hpp>

#include "SampleSets.hpp"
#include "cngpld/summarize.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

namespace {

struct SummaryRow {
	position pos;
	double value;
};

std::vector<SummaryRow> read_expected(const std::string& path) {
	std::ifstream in(path.c_str());
	BOOST_REQUIRE_MESSAGE(in.is_open(), "Failed to open " << path);
	std::string line;
	std::getline(in, line);
	std::vector<SummaryRow> rows;
	while (std::getline(in, line)) {
		if (line.empty()) continue;
		std::istringstream iss(line);
		SummaryRow row;
		iss >> row.pos >> row.value;
		rows.push_back(row);
	}
	return rows;
}

void check_summary(const cngpld::CNSummary& got, const std::vector<SummaryRow>& expected) {
	BOOST_REQUIRE_EQUAL(got.size(), expected.size());
	for (std::size_t i = 0; i < got.size(); ++i) {
		BOOST_CHECK_EQUAL(got[i].pos, expected[i].pos);
		BOOST_CHECK_CLOSE(got[i].value, expected[i].value, 1e-5);
	}
}

} // namespace

BOOST_AUTO_TEST_CASE(SummarizeCN_Matches_R_Case1_Amp_And_Del)
{
	cna::SegmentedSampleSet<rvalue> seg;
	seg.read("cngpld_case1_input.seg");
	check_summary(cngpld::summarize_cn(seg, "s1", 0, 1, 0.5), read_expected("cngpld_case1_amp_expected.tsv"));
	check_summary(cngpld::summarize_cn(seg, "s1", 0, -1, 0.5), read_expected("cngpld_case1_del_expected.tsv"));
}

BOOST_AUTO_TEST_CASE(SummarizeCN_Matches_R_Case2_Amp_And_Del)
{
	cna::SegmentedSampleSet<rvalue> seg;
	seg.read("cngpld_case2_input.seg");
	check_summary(cngpld::summarize_cn(seg, "s1", 0, 1, 0.5), read_expected("cngpld_case2_amp_expected.tsv"));
	check_summary(cngpld::summarize_cn(seg, "s1", 0, -1, 0.5), read_expected("cngpld_case2_del_expected.tsv"));
}

BOOST_AUTO_TEST_CASE(SummarizeCN_Preserves_Denominator_Semantics)
{
	cna::SegmentedSampleSet<rvalue> seg;
	seg.read("cngpld_case3_input.seg");
	check_summary(cngpld::summarize_cn(seg, "s1", 0, 1, 0.5), read_expected("cngpld_case3_amp_expected.tsv"));
}

BOOST_AUTO_TEST_CASE(SummarizeCN_Explicit_Positions_And_Empty_Overlap_Match_R)
{
	cna::SegmentedSampleSet<rvalue> seg;
	seg.read("cngpld_case4_input.seg");
	std::vector<position> positions;
	positions.push_back(50);
	positions.push_back(100);
	positions.push_back(120);
	positions.push_back(150);
	positions.push_back(200);
	positions.push_back(220);
	positions.push_back(250);
	check_summary(cngpld::summarize_cn(seg, "s1", 0, 1, 0.5, &positions), read_expected("cngpld_case4_amp_expected.tsv"));
}
