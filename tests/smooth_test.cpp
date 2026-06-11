#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Smooth CNA Tests"
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cbs/smooth.hpp"

using namespace std;

static ifstream open_test_file(const string& path) {
	ifstream in(path.c_str());
	if (!in.is_open()) in.open((string("tests/") + path).c_str());
	if (!in.is_open()) in.open((string("../tests/") + path).c_str());
	return in;
}

static vector<double> read_second_column_numeric(const string& path) {
	ifstream in = open_test_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<double> values;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string a, b;
		getline(ss, a, '\t');
		getline(ss, b, '\t');
		if (b == "NA" || b == "NaN") values.push_back(numeric_limits<double>::quiet_NaN());
		else if (b == "Inf") values.push_back(numeric_limits<double>::infinity());
		else if (b == "-Inf") values.push_back(-numeric_limits<double>::infinity());
		else values.push_back(stod(b));
	}
	return values;
}

static vector<int> read_second_column_int(const string& path) {
	ifstream in = open_test_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	vector<int> values;
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		string a, b;
		getline(ss, a, '\t');
		getline(ss, b, '\t');
		values.push_back(stoi(b));
	}
	return values;
}

static vector<vector<double>> read_matrix(const string& path) {
	ifstream in = open_test_file(path);
	BOOST_REQUIRE_MESSAGE(in.is_open(), string("Failed to open ") + path);
	string line;
	getline(in, line);
	istringstream hs(line);
	string field;
	vector<string> headers;
	while (getline(hs, field, '\t')) headers.push_back(field);
	vector<vector<double>> cols(headers.size() - 1);
	while (getline(in, line)) {
		if (line.empty()) continue;
		istringstream ss(line);
		getline(ss, field, '\t');
		for (size_t i = 1; i < headers.size(); ++i) {
			getline(ss, field, '\t');
			if (field == "NA" || field == "NaN") cols[i - 1].push_back(numeric_limits<double>::quiet_NaN());
			else if (field == "Inf") cols[i - 1].push_back(numeric_limits<double>::infinity());
			else if (field == "-Inf") cols[i - 1].push_back(-numeric_limits<double>::infinity());
			else cols[i - 1].push_back(stod(field));
		}
	}
	return cols;
}

BOOST_AUTO_TEST_SUITE(SmoothCNA)

BOOST_AUTO_TEST_CASE(FixtureFiles_AreReadable)
{
	auto x1 = read_second_column_numeric("smooth_case1_input.tsv");
	auto c1 = read_second_column_int("smooth_case1_chrom.tsv");
	auto y1 = read_second_column_numeric("smooth_case1_expected.tsv");
	BOOST_REQUIRE_EQUAL(x1.size(), c1.size());
	BOOST_REQUIRE_EQUAL(x1.size(), y1.size());
	BOOST_CHECK_EQUAL(c1.front(), 1);
	BOOST_CHECK_EQUAL(c1.back(), 2);
}

BOOST_AUTO_TEST_CASE(SingleProfile_Matches_DNAcopy)
{
	const auto x = read_second_column_numeric("smooth_case1_input.tsv");
	const auto chrom = read_second_column_int("smooth_case1_chrom.tsv");
	const auto expected = read_second_column_numeric("smooth_case1_expected.tsv");
	const auto observed = cbs::smooth(x, chrom);
	BOOST_REQUIRE_EQUAL(observed.size(), expected.size());
	for (size_t i = 0; i < observed.size(); ++i) {
		if (std::isfinite(expected[i])) BOOST_CHECK_SMALL(observed[i] - expected[i], 1e-9);
		else BOOST_CHECK(!std::isfinite(observed[i]));
	}
}

BOOST_AUTO_TEST_CASE(MissingValues_And_Boundaries_Match_DNAcopy)
{
	const auto x = read_second_column_numeric("smooth_case2_input.tsv");
	const auto chrom = read_second_column_int("smooth_case2_chrom.tsv");
	const auto expected = read_second_column_numeric("smooth_case2_expected.tsv");
	const auto observed = cbs::smooth(x, chrom, 2);
	BOOST_REQUIRE_EQUAL(observed.size(), expected.size());
	for (size_t i = 0; i < observed.size(); ++i) {
		if (std::isfinite(expected[i])) BOOST_CHECK_SMALL(observed[i] - expected[i], 1e-9);
		else BOOST_CHECK_EQUAL(std::isfinite(observed[i]), std::isfinite(expected[i]));
	}
}

BOOST_AUTO_TEST_CASE(MatrixSmoothing_Matches_DNAcopy)
{
	const auto samples = read_matrix("smooth_case3_input.tsv");
	const auto chrom = read_second_column_int("smooth_case3_chrom.tsv");
	const auto expected = read_matrix("smooth_case3_expected.tsv");
	const auto observed = cbs::smooth_matrix(samples, chrom, 2);
	BOOST_REQUIRE_EQUAL(observed.size(), expected.size());
	for (size_t j = 0; j < observed.size(); ++j) {
		BOOST_REQUIRE_EQUAL(observed[j].size(), expected[j].size());
		for (size_t i = 0; i < observed[j].size(); ++i) {
			if (std::isfinite(expected[j][i])) BOOST_CHECK_SMALL(observed[j][i] - expected[j][i], 1e-9);
			else BOOST_CHECK_EQUAL(std::isfinite(observed[j][i]), std::isfinite(expected[j][i]));
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
