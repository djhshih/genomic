#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Unit Tests for cna"
#include <boost/test/unit_test.hpp>

#include "config.h"
#include "global.hpp"
#include "SampleSets.hpp"
#include "FilesDiff.hpp"

#include <fstream>
#include <queue>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <sstream>

using namespace std;


class IOTests {
	friend struct IOFixture;
private:
	vector<string> configFilenames;
	IOTests() {
		configFilenames.push_back("iotest_basic.cfg");
		configFilenames.push_back("iotest_basic_as.cfg");
		configFilenames.push_back("iotest_conversion.cfg");
	}
};

struct IOFixture
{
	FilesDiff diff;
	fstream configf;
	IOTests tests;
	const char chrComment, chrClass;
	
	cna::RawSampleSet<rvalue> rset;
	cna::RawSampleSet<alleles_rcn> rset_as;
	cna::SegmentedSampleSet<rvalue> sset;
	cna::SegmentedSampleSet<alleles_rcn> sset_as;
	
	cna::GenericSampleSet gset;
	
	vector< queue<cna::SampleSet*> > sets;
	vector< queue<string> > filenames;
	
	void readConfigFile(string& filename, size_t type) {
		
		configf.open(filename.c_str(), ios::in);
		if (!configf.is_open()) {
			throw runtime_error("Failed to open IO test configuration file");
		}
		
		// fstream::eof() is not triggered until we read PAST the end of file
		// istringstream::eof() is triggered when we read TO the end of file (unless an extra return character is present)
		
		string line, s;
		istringstream stream;
		while (true) {
			getline(configf, line);
			if (configf.eof()) break;
			// ignore comments or blank lines
			if (line[0] == chrComment || line == "") continue;
			
			stream.clear();  // clear error flag
			stream.str(line);  // setup stream
		
			while (stream >> s) {
				if (s[0] == chrClass) {
					// item specifies which class to use
					// remove the first character
					s = s.substr(1);
					if (s == "RawSampleSet") {
						sets[type].push(&rset);
					} else if (s == "RawSampleSet<alleles_cn>" || s == "cna::RawSampleSet<alleles_cn>") {
						sets[type].push(&rset_as);
					} else if (s == "SegmentedSampleSet") {
						sets[type].push(&sset);
					} else if (s == "SegmentedSampleSet<alleles_cn>" || s == "cna::SegmentedSampleSet<alleles_cn>") {
						sets[type].push(&sset_as);
					} else if (s == "GenericSampleSet" || s == "cna::GenericSampleSet") {
						sets[type].push(&gset);
					} else {
						throw runtime_error("Invalid class specified");
					}
				} else {
					// item is a file name
					filenames[type].push(s);
				}
			}  // string stream loop
		}  // file loop
		
		configf.close();
	}
	
	void test() {
		for (size_t i = 0; i < sets.size(); ++i) {
			while (!sets[i].empty()) {
				cna::SampleSet* set = sets[i].front();
				sets[i].pop();
				set->read(filenames[i].front());
				filenames[i].pop();
				string out = filenames[i].front();
				filenames[i].pop();
				set->write(out);
				string ans = filenames[i].front();
				filenames[i].pop();
				BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
			}
			if (!filenames[i].empty()) {
				throw runtime_error("filenames queue was not exhausted: cna::SampleSet IO Test configuration file is likely malformed");
			}
		}
	}
	
	IOFixture() : diff(1), chrComment('#'), chrClass('@') {	
		
		sets.resize(tests.configFilenames.size());
		filenames.resize(tests.configFilenames.size());
		
		vector<string>::iterator it, end = tests.configFilenames.end();
		size_t type = -1;
		for (it = tests.configFilenames.begin(); it != end; ++it) {
			readConfigFile(*it, ++type);
		}
		
	}
		
	~IOFixture() {
		if (configf.is_open()) configf.close();
	}
};


namespace {

std::string shell_quote(const std::string& s) {
	std::string out = "'";
	for (char c : s) {
		if (c == '\'') out += "'\\''";
		else out += c;
	}
	out += "'";
	return out;
}

}

BOOST_AUTO_TEST_SUITE(SampleSetBasic)

BOOST_FIXTURE_TEST_CASE(InputOutput, IOFixture)
{
	BOOST_TEST_MESSAGE("Input output");
	test();
}

BOOST_AUTO_TEST_CASE(RawSampleSet_CopyConstructor)
{
	BOOST_TEST_MESSAGE("RawSampleSet copy constructor");
	
	FilesDiff diff;
	
	string fn_in = "raw1.in";
	string fn_copy = "raw1.copy";
	string fn_copy2 = "raw1.copy2";
	
	typedef cna::RawSampleSet<rvalue> SampleSetType;
	
	SampleSetType* set = new SampleSetType();
	set->read(fn_in);
	set->write(fn_copy);
	
	SampleSetType copy(*set);
	delete set;
	copy.write(fn_copy2);
	
	BOOST_CHECK_EQUAL(diff.different(fn_copy, fn_copy2), 0);
}

BOOST_AUTO_TEST_CASE(SegmentedSampleSet_CopyConstructor)
{
	BOOST_TEST_MESSAGE("SegmentedSampleSet copy constructor");
	
	FilesDiff diff;
	
	string fn_in = "seg1.in";
	string fn_copy = "seg1.copy";
	string fn_copy2 = "seg1.copy2";
	
	typedef cna::SegmentedSampleSet<rvalue> SampleSetType;
	
	SampleSetType* set = new SampleSetType();
	set->read(fn_in);
	set->write(fn_copy);
	
	SampleSetType copy(*set);
	delete set;
	copy.write(fn_copy2);
	
	BOOST_CHECK_EQUAL(diff.different(fn_copy, fn_copy2), 0);
}

BOOST_AUTO_TEST_CASE(SegmentedSampleSet_Filter)
{
	BOOST_TEST_MESSAGE("SegmentedSampleSet filter");
	
	FilesDiff diff;
	
	cna::SegmentedSampleSet<rvalue> set1a, set1b, set2a;
	set1a.read("segfilt1a.in");
	set2a.read("segfilt1a.in");
	set1b.read("segfilt1b.in");
	
	string out_seg = "segfilt1.out", ans_seg = "segfilt1.ans";
	string out2_seg = "segfilt2.out", ans2_seg = "segfilt2.ans";
	
	set1a.filter(set1b, 0.8);
	set1a.write(out_seg);
	
	// invert filter
	set2a.filter(set1b, 0.8, true);
	set2a.write(out2_seg);
	
	BOOST_CHECK_EQUAL(diff.different(out_seg, ans_seg), 0);
	BOOST_CHECK_EQUAL(diff.different(out2_seg, ans2_seg), 0);
}

BOOST_AUTO_TEST_CASE(RawSampleSet_Filter)
{
	BOOST_TEST_MESSAGE("RawSampleSet filter");
	
	FilesDiff diff;
	
	const unsigned ncases = 3;
	const char* set1fnames[] = {"cnfilt1a.in", "cnfilt2a.in", "cnfilt3a.in"};
	const char* set2fnames[] = {"cnfilt1b.in", "cnfilt2b.in", "cnfilt3b.in"};
	const char* outfnames[] = {"cnfilt1.out", "cnfilt2.out", "cnfilt3.out"};
	const char* ansfnames[] = {"cnfilt1.ans", "cnfilt2.ans", "cnfilt3.ans"};
	
	for (unsigned i = 0; i < ncases; i++) {
	
		cna::RawSampleSet<cnvalue> set1, set2;
		set1.read(set1fnames[i]);
		set2.read(set2fnames[i]);
		
		string out = outfnames[i], ans = ansfnames[i];
		
		set1.filter(set2);
		set1.write(out);
		
		BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
		
	}
}

BOOST_AUTO_TEST_CASE(CLI_Segment_Matches_DNAcopy_Expected_Output)
{
	FilesDiff diff;
	const std::string input = "segment_cli_case1_input.cn";
	const std::string output = "segment_cli_case1_output.seg";
	const std::string expected = "segment_cli_case1_expected.seg";

	const std::string cmd = std::string("../cna segment -i ") + shell_quote(input) +
		" -o " + shell_quote(output);
	const int rc = std::system(cmd.c_str());
	BOOST_REQUIRE_EQUAL(rc, 0);
	BOOST_CHECK_EQUAL(diff.different(output, expected), 0);
}

BOOST_AUTO_TEST_CASE(RawSampleSet_Filter_RemovesAlternatingMarkersAndPreservesAlignment)
{
	const char* set1fn = "cnfilt_cleanup_a.in";
	const char* set2fn = "cnfilt_cleanup_b.in";
	const char* outfn = "cnfilt_cleanup.out";
	{
		ofstream out(set1fn);
		out << "marker\tchromosome\tposition\ts1\ts2\n";
		out << "m1\tchr1\t10\t1\t10\n";
		out << "m2\tchr1\t20\t2\t20\n";
		out << "m3\tchr1\t30\t3\t30\n";
		out << "m4\tchr1\t40\t4\t40\n";
	}
	{
		ofstream out(set2fn);
		out << "marker\tchromosome\tposition\n";
		out << "m2\tchr1\t20\n";
		out << "m4\tchr1\t40\n";
	}

	cna::RawSampleSet<cnvalue> set1, set2;
	set1.read(string(set1fn));
	set2.read(string(set2fn));
	set1.filter(set2);
	set1.write(string(outfn));

	ifstream in(outfn);
	string line;
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "marker\tchromosome\tposition\ts1\ts2");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "m1\t1\t10\t1\t10");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "m3\t1\t30\t3\t30");
	BOOST_CHECK(!getline(in, line));

	std::remove(set1fn);
	std::remove(set2fn);
	std::remove(outfn);
}

BOOST_AUTO_TEST_CASE(SpuriousSegmentFilter_UsesStrictCountThreshold)
{
	cna::spurious_segment_filter<rvalue> filter(3);

	cna::Segment<rvalue> below(1, 1, 10, 2, 0.0f);
	cna::Segment<rvalue> atThreshold(1, 1, 10, 3, 0.0f);
	cna::Segment<rvalue> above(1, 1, 10, 4, 0.0f);

	BOOST_CHECK_EQUAL(filter(below), true);
	BOOST_CHECK_EQUAL(filter(atThreshold), false);
	BOOST_CHECK_EQUAL(filter(above), false);
}

BOOST_AUTO_TEST_CASE(SmallSegmentFilter_UsesStrictLengthThreshold)
{
	cna::small_segment_filter<rvalue> filter(10);

	cna::Segment<rvalue> below(1, 1, 9, 1, 0.0f);
	cna::Segment<rvalue> atThreshold(1, 1, 10, 1, 0.0f);
	cna::Segment<rvalue> above(1, 1, 11, 1, 0.0f);

	BOOST_CHECK_EQUAL(filter(below), true);
	BOOST_CHECK_EQUAL(filter(atThreshold), false);
	BOOST_CHECK_EQUAL(filter(above), false);
}

BOOST_AUTO_TEST_CASE(BalancedSegmentFilter_ScalarThresholds)
{
	cna::balanced_segment_filter<rvalue> filter(0.0f, 0.2f);

	cna::Segment<rvalue> inside(1, 1, 10, 1, 0.1f);
	cna::Segment<rvalue> lowerEdge(1, 1, 10, 1, -0.2f);
	cna::Segment<rvalue> upperEdge(1, 1, 10, 1, 0.2f);
	cna::Segment<rvalue> below(1, 1, 10, 1, -0.3f);
	cna::Segment<rvalue> above(1, 1, 10, 1, 0.3f);

	BOOST_CHECK_EQUAL(filter(inside), false);
	BOOST_CHECK_EQUAL(filter(lowerEdge), true);
	BOOST_CHECK_EQUAL(filter(upperEdge), true);
	BOOST_CHECK_EQUAL(filter(below), true);
	BOOST_CHECK_EQUAL(filter(above), true);
}

BOOST_AUTO_TEST_CASE(BalancedSegmentFilter_AlleleSpecificUsesTotalSignal)
{
	cna::balanced_segment_filter<alleles_cn> filter(2.0f, 0.3f);

	cna::Segment<alleles_cn> balanced(1, 1, 10, 1, alleles_cn(1, 1));
	cna::Segment<alleles_cn> highTotal(1, 1, 10, 1, alleles_cn(2, 1));
	cna::Segment<alleles_cn> lowTotal(1, 1, 10, 1, alleles_cn(1, 0));
	cna::Segment<alleles_cn> imbalancedButSameTotal(1, 1, 10, 1, alleles_cn(2, 0));

	BOOST_CHECK_EQUAL(filter(balanced), false);
	BOOST_CHECK_EQUAL(filter(highTotal), true);
	BOOST_CHECK_EQUAL(filter(lowTotal), true);
	BOOST_CHECK_EQUAL(filter(imbalancedButSameTotal), false);
}

BOOST_AUTO_TEST_CASE(BalancedSegmentFilter_AlleleSpecificRelativeCNUsesTotalSignal)
{
	cna::balanced_segment_filter<alleles_rcn> filter(0.0f, 0.2f);

	cna::Segment<alleles_rcn> balanced(1, 1, 10, 1, alleles_rcn(0, 0));
	cna::Segment<alleles_rcn> highTotal(1, 1, 10, 1, alleles_rcn(1, 0));
	cna::Segment<alleles_rcn> lowTotal(1, 1, 10, 1, alleles_rcn(-1, 0));
	cna::Segment<alleles_rcn> imbalancedButSameTotal(1, 1, 10, 1, alleles_rcn(1, -1));

	BOOST_CHECK_EQUAL(filter(balanced), false);
	BOOST_CHECK_EQUAL(filter(highTotal), true);
	BOOST_CHECK_EQUAL(filter(lowTotal), true);
	BOOST_CHECK_EQUAL(filter(imbalancedButSameTotal), false);
}

BOOST_AUTO_TEST_CASE(SegmentedSampleSet_Find_ExactAndLowerBound)
{
	cna::SegmentedSampleSet<rvalue> set;
	cna::SegmentedSampleSet<rvalue>::SegmentedSample* sample = set.create("sample1");

	cna::Segment<rvalue> seg1(1, 10, 19, 10, 1.0f);
	cna::Segment<rvalue> seg2(1, 20, 29, 10, 2.0f);
	cna::Segment<rvalue> seg3(1, 40, 49, 10, 3.0f);

	sample->addToChromosome(0, seg1);
	sample->addToChromosome(0, seg2);
	sample->addToChromosome(0, seg3);

	BOOST_CHECK_EQUAL(set.find(sample, 0, 10), 0u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 20), 1u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 40), 2u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 25), 1u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 39), 1u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 5), 0u);
	BOOST_CHECK_EQUAL(set.find(sample, 0, 100), 2u);
}

BOOST_AUTO_TEST_CASE(SegmentedSampleSet_Find_EmptyChromosome)
{
	cna::SegmentedSampleSet<rvalue> set;
	cna::SegmentedSampleSet<rvalue>::SegmentedSample* sample = set.create("sample1");

	BOOST_CHECK_EQUAL(set.find(sample, 0, 10), 0u);
}

BOOST_AUTO_TEST_CASE(RawSampleSet_Read_WithConfiguredDelimiter)
{
	const char* input = "raw_delim_test.in";
	const char* output = "raw_delim_test.out";
	{
		ofstream out(input);
		out << "marker,chromosome,position,s1,s2\n";
		out << "m1,chr1,10,1.5,2.5\n";
		out << "m2,chr1,20,3.5,4.5\n";
	}

	cna::RawSampleSet<rvalue> set;
	set.setIO(IOProperties(',', 1, 0, false));
	set.read(string(input));
	set.write(string(output));

	ifstream in(output);
	string line;
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "marker,chromosome,position,s1,s2");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "m1,1,10,1.5,2.5");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "m2,1,20,3.5,4.5");

	std::remove(input);
	std::remove(output);
}

BOOST_AUTO_TEST_CASE(SegmentedSampleSet_Read_WithConfiguredDelimiter)
{
	const char* input = "seg_delim_test.in";
	const char* output = "seg_delim_test.out";
	{
		ofstream out(input);
		out << "sample,chromosome,start,end,count,state\n";
		out << "s1,chr1,10,20,3,1.5\n";
		out << "s1,chr1,21,30,2,2.5\n";
	}

	cna::SegmentedSampleSet<rvalue> set;
	set.setIO(IOProperties(',', 1, 0, false));
	set.read(string(input));
	set.write(string(output));

	ifstream in(output);
	string line;
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "sample,chromosome,start,end,count,state");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "s1,1,10,20,3,1.5");
	getline(in, line);
	BOOST_CHECK_EQUAL(line, "s1,1,21,30,2,2.5");

	std::remove(input);
	std::remove(output);
}

BOOST_AUTO_TEST_CASE(RawSampleSet_InvalidInputPath)
{
	BOOST_CHECK_THROW(cna::RawSampleSet<rvalue>().read("does-not-exist.cn"), runtime_error);
}

BOOST_AUTO_TEST_CASE(GenericSampleSet_InvalidInputPath)
{
	BOOST_CHECK_THROW(cna::GenericSampleSet().read("does-not-exist.cn"), runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(Constants)

BOOST_AUTO_TEST_CASE(mapping)
{
	using namespace cna::mapping;
	BOOST_CHECK_EQUAL(chromosome["chr1"], 1);
	BOOST_CHECK_EQUAL(chromosome["chr9"], 9);
	BOOST_CHECK_EQUAL(chromosome["X"], 23);
	BOOST_CHECK_EQUAL(chromosome["chrY"], 24);
	BOOST_CHECK_EQUAL(chromosome["30"], 0);
}

BOOST_AUTO_TEST_SUITE_END()
