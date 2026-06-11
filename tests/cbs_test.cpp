#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "CBS Tests"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <tuple>

#include "cbs/CBS.hpp"

using namespace std;

struct ExpectedSegmentation {
	vector<int> lengths;
	vector<double> means;
};

static ifstream open_expected_file(const string& path) {
	ifstream in(path.c_str());
	if (!in.is_open()) in.open((string("tests/") + path).c_str());
	if (!in.is_open()) in.open((string("../tests/") + path).c_str());
	return in;
}

static vector<int> read_segment_starts(const string& path) {
	ifstream in = open_expected_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<int> starts;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		starts.push_back(std::stoi(field));
	}
	return starts;
}

static ExpectedSegmentation read_expected_segmentation(const string& path, int n) {
	ifstream in = open_expected_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<int> starts;
	vector<double> means;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		starts.push_back(std::stoi(field));
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		means.push_back(std::stod(field));
	}
	vector<int> lengths;
	for (size_t i = 0; i < starts.size(); ++i) {
		const int end = (i + 1 < starts.size()) ? (starts[i + 1] - 1) : n;
		lengths.push_back(end - starts[i] + 1);
	}
	return {lengths, means};
}

static vector<double> read_segment_means(const string& path) {
	ifstream in = open_expected_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<double> means;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		means.push_back(std::stod(field));
	}
	return means;
}

static vector<double> read_values_second_column(const string& path) {
	ifstream in = open_expected_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<double> values;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string field;
		getline(ss, field, '\t');
		getline(ss, field, '\t');
		values.push_back(std::stod(field));
	}
	return values;
}

BOOST_AUTO_TEST_SUITE(CBS)

BOOST_AUTO_TEST_CASE(DNAcopy_ExpectedFiles_ArePresentAndSane)
{
	vector<double> x1 = read_values_second_column("cbs_case1_input.tsv");
	BOOST_REQUIRE_EQUAL(x1.size(), 60u);
	BOOST_CHECK_SMALL(x1[0], 1e-12);
	BOOST_CHECK_CLOSE(x1[20], 1.5, 1e-9);
	BOOST_CHECK_SMALL(x1[40], 1e-12);

	vector<double> x2 = read_values_second_column("cbs_case2_weighted_input.tsv");
	vector<double> w2 = read_values_second_column("cbs_case2_weighted_weights.tsv");
	BOOST_REQUIRE_EQUAL(x2.size(), 60u);
	BOOST_REQUIRE_EQUAL(w2.size(), 60u);
	BOOST_CHECK_CLOSE(x2[15], 2.0, 1e-9);
	BOOST_CHECK_CLOSE(x2[30], -1.5, 1e-9);
	BOOST_CHECK_CLOSE(w2[15], 0.5, 1e-9);
	BOOST_CHECK_CLOSE(w2[30], 2.0, 1e-9);

	vector<double> x3 = read_values_second_column("cbs_case3_noisy_input.tsv");
	vector<double> x4 = read_values_second_column("cbs_case4_noisy_input.tsv");
	BOOST_REQUIRE_EQUAL(x3.size(), 100u);
	BOOST_REQUIRE_EQUAL(x4.size(), 100u);

	vector<int> starts1 = read_segment_starts("cbs_case1_expected.tsv");
	vector<double> means1 = read_segment_means("cbs_case1_expected.tsv");
	BOOST_REQUIRE_EQUAL(starts1.size(), 3u);
	BOOST_REQUIRE_EQUAL(means1.size(), 3u);
	BOOST_CHECK_EQUAL(starts1[0], 1);
	BOOST_CHECK_EQUAL(starts1[1], 21);
	BOOST_CHECK_EQUAL(starts1[2], 41);
	BOOST_CHECK_SMALL(means1[0], 1e-12);
	BOOST_CHECK_CLOSE(means1[1], 1.5, 1e-9);
	BOOST_CHECK_SMALL(means1[2], 1e-12);

	vector<int> starts2 = read_segment_starts("cbs_case2_weighted_expected.tsv");
	vector<double> means2 = read_segment_means("cbs_case2_weighted_expected.tsv");
	BOOST_REQUIRE_EQUAL(starts2.size(), 4u);
	BOOST_REQUIRE_EQUAL(means2.size(), 4u);
	BOOST_CHECK_EQUAL(starts2[0], 1);
	BOOST_CHECK_EQUAL(starts2[1], 16);
	BOOST_CHECK_EQUAL(starts2[2], 31);
	BOOST_CHECK_EQUAL(starts2[3], 46);
	BOOST_CHECK_SMALL(means2[0], 1e-12);
	BOOST_CHECK_CLOSE(means2[1], 2.0, 1e-9);
	BOOST_CHECK_CLOSE(means2[2], -1.5, 1e-9);
	BOOST_CHECK_SMALL(means2[3], 1e-12);
}

BOOST_AUTO_TEST_CASE(Unweighted_tmaxo_Matches_FortranRawStatisticBehavior)
{
	const vector<double> x = read_values_second_column("cbs_case1_input.tsv");
	double sumx = 0.0, sumsq = 0.0;
	for (double v : x) {
		sumx += v;
		sumsq += v * v;
	}
	const double tss = sumsq - (sumx * sumx) / x.size();
	const vector<tuple<int, bool, int, int, double>> cases{{2, false, 0, 58, 1000.0}, {2, true, 0, 58, 10.0}};
	for (const auto& tc : cases) {
		const int al0 = std::get<0>(tc);
		const bool ibin = std::get<1>(tc);
		const int expect_start = std::get<2>(tc);
		const int expect_end = std::get<3>(tc);
		const double min_stat = std::get<4>(tc);
		const auto obs = cbs::tmaxo(x, tss, al0, ibin);
		BOOST_TEST_CONTEXT("al0=" << al0 << " ibin=" << ibin) {
			BOOST_CHECK_EQUAL(obs.start, expect_start);
			BOOST_CHECK_EQUAL(obs.end, expect_end);
			BOOST_CHECK(obs.statistic > min_stat);
		}
	}
}

BOOST_AUTO_TEST_CASE(Weighted_wtmaxo_Matches_FortranRawStatisticBehavior)
{
	const vector<double> x = read_values_second_column("cbs_case2_weighted_input.tsv");
	const vector<double> wts = read_values_second_column("cbs_case2_weighted_weights.tsv");
	vector<double> cwts(wts.size());
	double wsum = 0.0, wxsum = 0.0, wxxsum = 0.0, csum = 0.0;
	for (size_t i = 0; i < x.size(); ++i) {
		csum += wts[i];
		cwts[i] = csum;
		wsum += wts[i];
		wxsum += wts[i] * x[i];
		wxxsum += wts[i] * x[i] * x[i];
	}
	const double tss = wxxsum - wsum * std::pow(wxsum / wsum, 2.0);
	const vector<tuple<int, int, int>> cases{{2, 0, 58}, {3, 0, 57}};
	for (const auto& tc : cases) {
		const int al0 = std::get<0>(tc);
		const auto obs = cbs::wtmaxo(x, wts, tss, cwts, al0);
		BOOST_TEST_CONTEXT("al0=" << al0) {
			BOOST_CHECK_EQUAL(obs.start, std::get<1>(tc));
			BOOST_CHECK_EQUAL(obs.end, std::get<2>(tc));
			BOOST_CHECK(obs.statistic > 1.0);
		}
	}
}

BOOST_AUTO_TEST_CASE(NoisyProfiles_ExpectedFiles_AreReadable)
{
	const vector<double> x3 = read_values_second_column("cbs_case3_noisy_input.tsv");
	const vector<double> x4 = read_values_second_column("cbs_case4_noisy_input.tsv");
	const vector<int> starts3 = read_segment_starts("cbs_case3_noisy_expected.tsv");
	const vector<int> starts4 = read_segment_starts("cbs_case4_noisy_expected.tsv");
	const vector<double> means3 = read_segment_means("cbs_case3_noisy_expected.tsv");
	const vector<double> means4 = read_segment_means("cbs_case4_noisy_expected.tsv");

	BOOST_REQUIRE_EQUAL(x3.size(), 100u);
	BOOST_REQUIRE_EQUAL(x4.size(), 100u);
	BOOST_CHECK(!starts3.empty());
	BOOST_CHECK(!starts4.empty());
	BOOST_REQUIRE_EQUAL(starts3.size(), means3.size());
	BOOST_REQUIRE_EQUAL(starts4.size(), means4.size());
}

BOOST_AUTO_TEST_CASE(Unweighted_SegmentDriver_MatchesDNAcopy_NoisyCase3)
{
	const vector<double> x = read_values_second_column("cbs_case3_noisy_input.tsv");
	const vector<tuple<string, double, int, bool, int>> cases{
		{"cbs_case3_noisy_expected.tsv", 0.01, 200, false, 2},
		{"cbs_case3_noisy_perm_alt_expected.tsv", 0.05, 100, false, 3}
	};
	for (const auto& tc : cases) {
		const auto expected = read_expected_segmentation(std::get<0>(tc), static_cast<int>(x.size()));
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((std::get<1>(tc) > 0 ? std::get<2>(tc) + 1 : 201) * (std::get<1>(tc) > 0 ? std::get<2>(tc) + 2 : 202) / 2 + 2, std::get<2>(tc) + 1);
		const auto seg = cbs::segment(x, false, std::get<1>(tc), std::get<2>(tc), std::get<3>(tc), std::get<4>(tc), 25, 200, 0.05, sbdry, 1e-6, rng, false, 0.05);
		BOOST_TEST_CONTEXT("expected=" << std::get<0>(tc) << " alpha=" << std::get<1>(tc) << " nperm=" << std::get<2>(tc) << " hybrid=" << std::get<3>(tc) << " min_width=" << std::get<4>(tc)) {
			BOOST_REQUIRE_EQUAL(seg.lengths.size(), expected.lengths.size());
			BOOST_REQUIRE_EQUAL(seg.means.size(), expected.means.size());
			BOOST_CHECK_EQUAL_COLLECTIONS(seg.lengths.begin(), seg.lengths.end(), expected.lengths.begin(), expected.lengths.end());
			for (size_t i = 0; i < expected.means.size(); ++i) BOOST_CHECK_SMALL(seg.means[i] - expected.means[i], 5e-5);
		}
	}
}

BOOST_AUTO_TEST_CASE(Unweighted_SegmentDriver_MatchesDNAcopy_NoisyCase4)
{
	const vector<double> x = read_values_second_column("cbs_case4_noisy_input.tsv");
	const vector<tuple<string, double, int, bool, int>> cases{
		{"cbs_case4_noisy_expected.tsv", 0.01, 200, false, 2},
		{"cbs_case4_noisy_perm_alt_expected.tsv", 0.05, 100, false, 3}
	};
	for (const auto& tc : cases) {
		const auto expected = read_expected_segmentation(std::get<0>(tc), static_cast<int>(x.size()));
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((std::get<2>(tc) + 1) * (std::get<2>(tc) + 2) / 2 + 2, std::get<2>(tc) + 1);
		const auto seg = cbs::segment(x, false, std::get<1>(tc), std::get<2>(tc), std::get<3>(tc), std::get<4>(tc), 25, 200, 0.05, sbdry, 1e-6, rng, false, 0.05);
		BOOST_TEST_CONTEXT("expected=" << std::get<0>(tc) << " alpha=" << std::get<1>(tc) << " nperm=" << std::get<2>(tc) << " hybrid=" << std::get<3>(tc) << " min_width=" << std::get<4>(tc)) {
			BOOST_REQUIRE_EQUAL(seg.lengths.size(), expected.lengths.size());
			BOOST_REQUIRE_EQUAL(seg.means.size(), expected.means.size());
			BOOST_CHECK_EQUAL_COLLECTIONS(seg.lengths.begin(), seg.lengths.end(), expected.lengths.begin(), expected.lengths.end());
			for (size_t i = 0; i < expected.means.size(); ++i) BOOST_CHECK_SMALL(seg.means[i] - expected.means[i], 5e-5);
		}
	}
}

BOOST_AUTO_TEST_CASE(Unweighted_SegmentDriver_MatchesDNAcopy_NoisyHybrid)
{
	const vector<tuple<string, string, double, int, int>> cases{
		{"cbs_case3_noisy_input.tsv", "cbs_case3_noisy_hybrid_expected.tsv", 0.01, 200, 2},
		{"cbs_case3_noisy_input.tsv", "cbs_case3_noisy_hybrid_alt_expected.tsv", 0.05, 100, 3},
		{"cbs_case4_noisy_input.tsv", "cbs_case4_noisy_hybrid_expected.tsv", 0.01, 200, 2},
		{"cbs_case4_noisy_input.tsv", "cbs_case4_noisy_hybrid_alt_expected.tsv", 0.05, 100, 3}
	};
	for (const auto& tc : cases) {
		const vector<double> x = read_values_second_column(std::get<0>(tc));
		const auto expected = read_expected_segmentation(std::get<1>(tc), static_cast<int>(x.size()));
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((std::get<3>(tc) + 1) * (std::get<3>(tc) + 2) / 2 + 2, std::get<3>(tc) + 1);
		const auto seg = cbs::segment(x, false, std::get<2>(tc), std::get<3>(tc), true, std::get<4>(tc), 25, 200, 0.05, sbdry, 1e-6, rng, false, 0.05);
		BOOST_TEST_CONTEXT("input=" << std::get<0>(tc) << " expected=" << std::get<1>(tc) << " alpha=" << std::get<2>(tc) << " nperm=" << std::get<3>(tc) << " hybrid=true min_width=" << std::get<4>(tc)) {
			BOOST_REQUIRE_EQUAL(seg.lengths.size(), expected.lengths.size());
			BOOST_REQUIRE_EQUAL(seg.means.size(), expected.means.size());
			BOOST_CHECK_EQUAL_COLLECTIONS(seg.lengths.begin(), seg.lengths.end(), expected.lengths.begin(), expected.lengths.end());
			for (size_t i = 0; i < expected.means.size(); ++i) BOOST_CHECK_SMALL(seg.means[i] - expected.means[i], 5e-5);
		}
	}
}

BOOST_AUTO_TEST_CASE(Unweighted_SegmentDriver_MatchesDNAcopy_SimpleCase)
{
	const vector<double> x = read_values_second_column("cbs_case1_input.tsv");
	const vector<tuple<string, double, int, bool, int>> cases{
		{"cbs_case1_expected.tsv", 0.01, 200, false, 2},
		{"cbs_case1_perm_alt_expected.tsv", 0.05, 100, false, 3},
		{"cbs_case1_hybrid_expected.tsv", 0.01, 200, true, 2},
		{"cbs_case1_hybrid_alt_expected.tsv", 0.05, 100, true, 3}
	};
	for (const auto& tc : cases) {
		const auto expected = read_expected_segmentation(std::get<0>(tc), static_cast<int>(x.size()));
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((std::get<2>(tc) + 1) * (std::get<2>(tc) + 2) / 2 + 2, std::get<2>(tc) + 1);
		const auto seg = cbs::segment(x, false, std::get<1>(tc), std::get<2>(tc), std::get<3>(tc), std::get<4>(tc), 25, 200, 0.05, sbdry, 1e-6, rng, false, 0.05);
		BOOST_TEST_CONTEXT("expected=" << std::get<0>(tc) << " alpha=" << std::get<1>(tc) << " nperm=" << std::get<2>(tc) << " hybrid=" << std::get<3>(tc) << " min_width=" << std::get<4>(tc)) {
			BOOST_REQUIRE_EQUAL(seg.lengths.size(), expected.lengths.size());
			BOOST_REQUIRE_EQUAL(seg.means.size(), expected.means.size());
			BOOST_CHECK_EQUAL_COLLECTIONS(seg.lengths.begin(), seg.lengths.end(), expected.lengths.begin(), expected.lengths.end());
			for (size_t i = 0; i < expected.means.size(); ++i) BOOST_CHECK_SMALL(seg.means[i] - expected.means[i], 1e-9);
		}
	}
}

BOOST_AUTO_TEST_CASE(Weighted_SegmentDriver_MatchesDNAcopy_SimpleCase)
{
	const vector<double> x = read_values_second_column("cbs_case2_weighted_input.tsv");
	const vector<double> wts = read_values_second_column("cbs_case2_weighted_weights.tsv");
	const vector<tuple<string, double, int, bool, int>> cases{
		{"cbs_case2_weighted_expected.tsv", 0.01, 200, false, 2},
		{"cbs_case2_weighted_perm_alt_expected.tsv", 0.05, 100, false, 3},
		{"cbs_case2_weighted_hybrid_expected.tsv", 0.01, 200, true, 2},
		{"cbs_case2_weighted_hybrid_alt_expected.tsv", 0.05, 100, true, 3}
	};
	for (const auto& tc : cases) {
		const auto expected = read_expected_segmentation(std::get<0>(tc), static_cast<int>(x.size()));
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((std::get<2>(tc) + 1) * (std::get<2>(tc) + 2) / 2 + 2, std::get<2>(tc) + 1);
		const auto seg = cbs::segment_weighted(x, wts, std::get<1>(tc), std::get<2>(tc), std::get<3>(tc), std::get<4>(tc), 25, 200, 0.05, sbdry, 1e-6, rng, false, 0.05);
		BOOST_TEST_CONTEXT("expected=" << std::get<0>(tc) << " alpha=" << std::get<1>(tc) << " nperm=" << std::get<2>(tc) << " hybrid=" << std::get<3>(tc) << " min_width=" << std::get<4>(tc)) {
			BOOST_REQUIRE_EQUAL(seg.lengths.size(), expected.lengths.size());
			BOOST_REQUIRE_EQUAL(seg.means.size(), expected.means.size());
			BOOST_CHECK_EQUAL_COLLECTIONS(seg.lengths.begin(), seg.lengths.end(), expected.lengths.begin(), expected.lengths.end());
			for (size_t i = 0; i < expected.means.size(); ++i) BOOST_CHECK_SMALL(seg.means[i] - expected.means[i], 1e-9);
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()
