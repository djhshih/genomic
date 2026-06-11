#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "NCList Tests"
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "NCList.hpp"

using namespace std;

static ifstream open_nclist_file(const string& path) {
	ifstream in(path.c_str());
	if (!in.is_open()) in.open((string("tests/") + path).c_str());
	return in;
}

static vector<cna::NCList::IntervalRef> read_intervals(const string& path) {
	ifstream in = open_nclist_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<cna::NCList::IntervalRef> intervals;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		size_t id = static_cast<size_t>(std::stoul(field));
		getline(ss, field, '\t');
		position start = static_cast<position>(std::stoul(field));
		getline(ss, field, '\t');
		position end = static_cast<position>(std::stoul(field));
		intervals.push_back(cna::NCList::IntervalRef{start, end, id});
	}
	return intervals;
}

struct QueryRange {
	position start;
	position end;
};

static vector<QueryRange> read_queries(const string& path) {
	ifstream in = open_nclist_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<QueryRange> queries;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		position start = static_cast<position>(std::stoul(field));
		getline(ss, field, '\t');
		position end = static_cast<position>(std::stoul(field));
		queries.push_back(QueryRange{start, end});
	}
	return queries;
}

static vector<vector<size_t> > read_expected_hits(const string& path, size_t nqueries) {
	ifstream in = open_nclist_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<vector<size_t> > hits(nqueries);
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		size_t qid = static_cast<size_t>(std::stoul(field));
		getline(ss, field, '\t');
		size_t hit = static_cast<size_t>(std::stoul(field));
		hits[qid].push_back(hit);
	}
	return hits;
}

static void check_iranges_case(const string& stem) {
	const vector<cna::NCList::IntervalRef> intervals = read_intervals(stem + "_intervals.tsv");
	const vector<QueryRange> queries = read_queries(stem + "_queries.tsv");
	const vector<vector<size_t> > expected = read_expected_hits(stem + "_expected.tsv", queries.size());
	cna::NCList index(intervals);

	for (size_t qi = 0; qi < queries.size(); ++qi) {
		vector<size_t> hits;
		index.findOverlaps(queries[qi].start, queries[qi].end, back_inserter(hits));
		sort(hits.begin(), hits.end());
		vector<size_t> sorted_expected = expected[qi];
		sort(sorted_expected.begin(), sorted_expected.end());
		BOOST_CHECK_EQUAL(index.overlapsAny(queries[qi].start, queries[qi].end), !expected[qi].empty());
		BOOST_CHECK(hits == sorted_expected);
	}
}

static void check_iranges_case_order(const string& stem) {
	const vector<cna::NCList::IntervalRef> intervals = read_intervals(stem + "_intervals.tsv");
	const vector<QueryRange> queries = read_queries(stem + "_queries.tsv");
	const vector<vector<size_t> > expected = read_expected_hits(stem + "_expected.tsv", queries.size());
	cna::NCList index(intervals);

	for (size_t qi = 0; qi < queries.size(); ++qi) {
		vector<size_t> hits;
		index.findOverlaps(queries[qi].start, queries[qi].end, back_inserter(hits));
		BOOST_CHECK(hits == expected[qi]);
	}
}

BOOST_AUTO_TEST_SUITE(NCList)

BOOST_AUTO_TEST_CASE(EmptyIndex)
{
	cna::NCList index;
	std::vector<size_t> hits;
	index.findOverlaps(10, 20, std::back_inserter(hits));
	BOOST_CHECK(hits.empty());
	BOOST_CHECK_EQUAL(index.overlapsAny(10, 20), false);
}

BOOST_AUTO_TEST_CASE(BoundaryTouchAndDuplicates)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{10, 20, 7});
	intervals.push_back(cna::NCList::IntervalRef{20, 30, 8});
	intervals.push_back(cna::NCList::IntervalRef{20, 30, 9});
	cna::NCList index(intervals);

	std::vector<size_t> hits;
	index.findOverlaps(20, 20, std::back_inserter(hits));
	std::sort(hits.begin(), hits.end());
	BOOST_REQUIRE_EQUAL(hits.size(), 3u);
	BOOST_CHECK_EQUAL(hits[0], 7u);
	BOOST_CHECK_EQUAL(hits[1], 8u);
	BOOST_CHECK_EQUAL(hits[2], 9u);
	BOOST_CHECK_EQUAL(index.overlapsAny(31, 40), false);
	BOOST_CHECK_EQUAL(index.overlapsAny(30, 30), true);
}

BOOST_AUTO_TEST_CASE(NestedIntervals)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{1, 100, 0});
	intervals.push_back(cna::NCList::IntervalRef{10, 90, 1});
	intervals.push_back(cna::NCList::IntervalRef{20, 80, 2});
	intervals.push_back(cna::NCList::IntervalRef{200, 250, 3});
	cna::NCList index(intervals);

	std::vector<size_t> hits;
	index.findOverlaps(50, 60, std::back_inserter(hits));
	std::sort(hits.begin(), hits.end());
	BOOST_REQUIRE_EQUAL(hits.size(), 3u);
	BOOST_CHECK_EQUAL(hits[0], 0u);
	BOOST_CHECK_EQUAL(hits[1], 1u);
	BOOST_CHECK_EQUAL(hits[2], 2u);
}

BOOST_AUTO_TEST_CASE(ZeroLengthInterval)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{5, 5, 42});
	cna::NCList index(intervals);
	BOOST_CHECK_EQUAL(index.overlapsAny(5, 5), true);
	BOOST_CHECK_EQUAL(index.overlapsAny(6, 6), false);
}

BOOST_AUTO_TEST_CASE(InvalidBuildThrows)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{10, 9, 0});
	BOOST_CHECK_THROW(cna::NCList index(intervals), std::logic_error);
}

BOOST_AUTO_TEST_CASE(InvalidQueryThrows)
{
	cna::NCList index;
	BOOST_CHECK_THROW(index.overlapsAny(10, 9), std::logic_error);
	std::vector<size_t> hits;
	BOOST_CHECK_THROW(index.findOverlaps(10, 9, std::back_inserter(hits)), std::logic_error);
}

BOOST_AUTO_TEST_CASE(LargeDuplicateRanges)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	for (size_t i = 0; i < 50; ++i)
		intervals.push_back(cna::NCList::IntervalRef{100, 200, i});
	cna::NCList index(intervals);
	std::vector<size_t> hits;
	index.findOverlaps(150, 150, std::back_inserter(hits));
	BOOST_REQUIRE_EQUAL(hits.size(), 50u);
	std::sort(hits.begin(), hits.end());
	for (size_t i = 0; i < 50; ++i)
		BOOST_CHECK_EQUAL(hits[i], i);
}

BOOST_AUTO_TEST_CASE(InvariantsHold)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{1, 100, 0});
	intervals.push_back(cna::NCList::IntervalRef{10, 20, 1});
	intervals.push_back(cna::NCList::IntervalRef{30, 40, 2});
	intervals.push_back(cna::NCList::IntervalRef{35, 38, 3});
	intervals.push_back(cna::NCList::IntervalRef{200, 210, 4});
	cna::NCList index(intervals);
	cna::NCList::InvariantSummary summary = index.checkInvariants();
	BOOST_CHECK_EQUAL(summary.sibling_starts_monotone, true);
	BOOST_CHECK_EQUAL(summary.sibling_ends_monotone, true);
	BOOST_CHECK_EQUAL(summary.parent_contains_children, true);
}

BOOST_AUTO_TEST_CASE(InvariantsHoldWithSameStart)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{10, 100, 0});
	intervals.push_back(cna::NCList::IntervalRef{10, 90, 1});
	intervals.push_back(cna::NCList::IntervalRef{10, 80, 2});
	intervals.push_back(cna::NCList::IntervalRef{20, 30, 3});
	cna::NCList index(intervals);
	cna::NCList::InvariantSummary summary = index.checkInvariants();
	BOOST_CHECK_EQUAL(summary.sibling_starts_monotone, true);
	BOOST_CHECK_EQUAL(summary.sibling_ends_monotone, true);
	BOOST_CHECK_EQUAL(summary.parent_contains_children, true);
}

BOOST_AUTO_TEST_CASE(EmulatesIRangesTraversalOrderAndPruning)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{1, 10, 0});
	intervals.push_back(cna::NCList::IntervalRef{2, 3, 1});
	intervals.push_back(cna::NCList::IntervalRef{4, 5, 2});
	intervals.push_back(cna::NCList::IntervalRef{20, 30, 3});
	intervals.push_back(cna::NCList::IntervalRef{21, 22, 4});
	cna::NCList index(intervals);

	std::vector<size_t> hits;
	index.findOverlaps(4, 4, std::back_inserter(hits));
	BOOST_REQUIRE_EQUAL(hits.size(), 2u);
	BOOST_CHECK_EQUAL(hits[0], 0u);
	BOOST_CHECK_EQUAL(hits[1], 2u);

	hits.clear();
	index.findOverlaps(23, 23, std::back_inserter(hits));
	BOOST_REQUIRE_EQUAL(hits.size(), 1u);
	BOOST_CHECK_EQUAL(hits[0], 3u);
}

BOOST_AUTO_TEST_CASE(MatchesNaiveClosedIntervalOverlap)
{
	std::vector<cna::NCList::IntervalRef> intervals;
	intervals.push_back(cna::NCList::IntervalRef{1, 100, 0});
	intervals.push_back(cna::NCList::IntervalRef{10, 20, 1});
	intervals.push_back(cna::NCList::IntervalRef{15, 25, 2});
	intervals.push_back(cna::NCList::IntervalRef{30, 30, 3});
	intervals.push_back(cna::NCList::IntervalRef{50, 60, 4});
	intervals.push_back(cna::NCList::IntervalRef{70, 90, 5});
	cna::NCList index(intervals);

	for (position qstart = 1; qstart <= 95; qstart += 7) {
		for (position qend = qstart; qend <= 100; qend += 11) {
			std::vector<size_t> hits;
			index.findOverlaps(qstart, qend, std::back_inserter(hits));
			std::sort(hits.begin(), hits.end());

			std::vector<size_t> expected;
			for (size_t i = 0; i < intervals.size(); ++i) {
				if (!(intervals[i].end < qstart || intervals[i].start > qend))
					expected.push_back(intervals[i].id);
			}
			std::sort(expected.begin(), expected.end());
			BOOST_CHECK(hits == expected);
			BOOST_CHECK_EQUAL(index.overlapsAny(qstart, qend), !expected.empty());
		}
	}
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase1)
{
	check_iranges_case("nclist_case1");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesOrderCase1)
{
	check_iranges_case_order("nclist_case1");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesOrderCase2)
{
	check_iranges_case_order("nclist_case2");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase2)
{
	check_iranges_case("nclist_case2");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase3)
{
	check_iranges_case("nclist_case3");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase4Random)
{
	check_iranges_case("nclist_case4");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase10Duplicates)
{
	check_iranges_case("nclist_case10_duplicates");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase5DeepContainment)
{
	check_iranges_case("nclist_case5_deep");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase6SameStart)
{
	check_iranges_case("nclist_case6_same_start");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase7SameEnd)
{
	check_iranges_case("nclist_case7_same_end");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase8Sparse)
{
	check_iranges_case("nclist_case8_sparse");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesCase9Dense)
{
	check_iranges_case("nclist_case9_dense");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesFuzz1)
{
	check_iranges_case("nclist_fuzz1");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesFuzz2)
{
	check_iranges_case("nclist_fuzz2");
}

BOOST_AUTO_TEST_CASE(MatchesIRangesFuzz3)
{
	check_iranges_case("nclist_fuzz3");
}

BOOST_AUTO_TEST_SUITE_END()
