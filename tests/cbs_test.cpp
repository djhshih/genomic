#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "CBS Tests"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>

#include "cbs/CBS.hpp"

using namespace std;

static ifstream open_expected_file(const string& path) {
	ifstream in(path.c_str());
	if (!in.is_open()) in.open((string("tests/") + path).c_str());
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

BOOST_AUTO_TEST_CASE(Unweighted_PortMatchesDNAcopy)
{
	const vector<double> x = read_values_second_column("cbs_case1_input.tsv");
	const vector<int> starts = read_segment_starts("cbs_case1_expected.tsv");
	const vector<double> means = read_segment_means("cbs_case1_expected.tsv");

	double sumx = 0.0, sumsq = 0.0;
	for (double v : x) {
		sumx += v;
		sumsq += v * v;
	}
	const double tss = sumsq - (sumx * sumx) / x.size();
	const auto obs = cbs::tmaxo(x, tss, 2, false);

	BOOST_REQUIRE_EQUAL(starts.size(), 3u);
	BOOST_REQUIRE_EQUAL(means.size(), 3u);
	BOOST_CHECK_EQUAL(obs.start, starts[1] - 1);
	BOOST_CHECK_EQUAL(obs.end, starts[2] - 1);
	BOOST_CHECK_CLOSE(means[1], 1.5, 1e-9);
}

BOOST_AUTO_TEST_CASE(Weighted_PortMatchesDNAcopy)
{
	const vector<double> x = read_values_second_column("cbs_case2_weighted_input.tsv");
	const vector<double> wts = read_values_second_column("cbs_case2_weighted_weights.tsv");
	const vector<int> starts = read_segment_starts("cbs_case2_weighted_expected.tsv");
	const vector<double> means = read_segment_means("cbs_case2_weighted_expected.tsv");

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
	const auto obs = cbs::wtmaxo(x, wts, tss, cwts, 2);

	BOOST_REQUIRE_EQUAL(starts.size(), 4u);
	BOOST_REQUIRE_EQUAL(means.size(), 4u);
	BOOST_CHECK_EQUAL(obs.start, starts[1] - 1);
	BOOST_CHECK_EQUAL(obs.end, starts[3] - 1);
	BOOST_CHECK_CLOSE(means[1], 2.0, 1e-9);
	BOOST_CHECK_CLOSE(means[2], -1.5, 1e-9);
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

BOOST_AUTO_TEST_CASE(NoisyProfile1_PortBoundaryProbeMatchesDNAcopy)
{
	const vector<double> x = read_values_second_column("cbs_case3_noisy_input.tsv");
	const vector<int> starts = read_segment_starts("cbs_case3_noisy_expected.tsv");
	const vector<double> means = read_segment_means("cbs_case3_noisy_expected.tsv");

	double sumx = 0.0, sumsq = 0.0;
	for (double v : x) {
		sumx += v;
		sumsq += v * v;
	}
	const double tss = sumsq - (sumx * sumx) / x.size();
	const auto obs = cbs::tmaxo(x, tss, 2, false);

	BOOST_REQUIRE_EQUAL(starts.size(), means.size());
	BOOST_REQUIRE(starts.size() >= 2u);
	BOOST_CHECK_EQUAL(obs.start, starts[1] - 1);
	BOOST_CHECK(obs.end >= obs.start);
	BOOST_CHECK(obs.end < static_cast<int>(x.size()));
}

BOOST_AUTO_TEST_CASE(NoisyProfile2_PortBoundaryProbeMatchesDNAcopy)
{
	const vector<double> x = read_values_second_column("cbs_case4_noisy_input.tsv");
	const vector<int> starts = read_segment_starts("cbs_case4_noisy_expected.tsv");
	const vector<double> means = read_segment_means("cbs_case4_noisy_expected.tsv");

	double sumx = 0.0, sumsq = 0.0;
	for (double v : x) {
		sumx += v;
		sumsq += v * v;
	}
	const double tss = sumsq - (sumx * sumx) / x.size();
	const auto obs = cbs::tmaxo(x, tss, 2, false);

	BOOST_REQUIRE_EQUAL(starts.size(), means.size());
	BOOST_REQUIRE(starts.size() >= 2u);
	BOOST_CHECK_EQUAL(obs.start, starts[1] - 1);
	BOOST_CHECK(obs.end >= obs.start);
	BOOST_CHECK(obs.end < static_cast<int>(x.size()));
}

BOOST_AUTO_TEST_SUITE_END()
