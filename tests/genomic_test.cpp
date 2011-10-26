#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Unit Tests for genomic"
#include <boost/test/unit_test.hpp>

#include "config.h"
#include "../global.hpp"
#include "../SampleSets.hpp"
#include "FilesDiff.hpp"

#include <fstream>
#include <queue>
#include <stdexcept>

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
	
	RawSampleSet<rvalue> rset;
	RawSampleSet<alleles_rcn> rset_as;
	SegmentedSampleSet<rvalue> sset;
	SegmentedSampleSet<alleles_rcn> sset_as;
	
	GenericSampleSet gset;
	
	vector< queue<SampleSet*> > sets;
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
		int count = 0;
		while (true) {
			getline(configf, line);
			if (configf.eof()) break;
			// ignore comments or blank lines
			if (line[0] == chrComment || line == "") continue;
			
			stream.clear();  // clear error flag
			stream.str(line);  // setup stream
		
			while (!stream.eof()) {
				stream >> s;
				if (s[0] == chrClass) {
					// item specifies which class to use
					// remove the first character
					s = s.substr(1);
					if (s == "RawSampleSet") {
						sets[type].push(&rset);
					} else if (s == "RawSampleSet<alleles_cn>") {
						sets[type].push(&rset_as);
					} else if (s == "SegmentedSampleSet") {
						sets[type].push(&sset);
					} else if (s == "SegmentedSampleSet<alleles_cn>") {
						sets[type].push(&sset_as);
					} else if (s == "GenericSampleSet") {
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
				SampleSet* set = sets[i].front();
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
				throw runtime_error("filenames queue was not exhausted: SampleSet IO Test configuration file is likely malformed");
			}
		}
	}
	
	IOFixture() : chrComment('#'), chrClass('@'), diff(1) {	
		
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


BOOST_AUTO_TEST_SUITE(SampleSetBasic)

BOOST_FIXTURE_TEST_CASE(InputOutput, IOFixture)
{
	BOOST_TEST_MESSAGE("Input output");
	test();
}

BOOST_AUTO_TEST_CASE(InputOutput_Picnic)
{
	BOOST_TEST_MESSAGE("Input output - Picnic");
	
	FilesDiff diff;
	
	vector<string> fileNames(2);
	fileNames[0] = "picnic1a.in";
	fileNames[1] = "picnic1b.in";
	string markersFileName = "picnic.snp6.in";
	string out = "picnic1.out";
	string ans = "picnic1.ans";
	string out_seg = "picnic1.seg";
	string ans_seg = "picnic1.seg.ans";
	
	PicnicSampleSet pset;
	pset.read(fileNames, markersFileName);
	pset.write(out);
	BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
	
	SegmentedSampleSet<PicnicSampleSet::Value> sset(pset);
	sset.write(out_seg);
	BOOST_CHECK_EQUAL(diff.different(out_seg, ans_seg), 0);
}

BOOST_AUTO_TEST_CASE(RawSampleSet_CopyConstructor)
{
	BOOST_TEST_MESSAGE("RawSampleSet copy constructor");
	
	FilesDiff diff;
	
	string fn_in = "raw1.in";
	string fn_copy = "raw1.copy";
	string fn_copy2 = "raw1.copy2";
	
	typedef RawSampleSet<rvalue> SampleSetType;
	
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
	
	typedef SegmentedSampleSet<rvalue> SampleSetType;
	
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
	
	SegmentedSampleSet<rvalue> set1a, set1b, set2a;
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
	
		RawSampleSet<cnvalue> set1, set2;
		set1.read(set1fnames[i]);
		set2.read(set2fnames[i]);
		
		string out = outfnames[i], ans = ansfnames[i];
		
		set1.filter(set2);
		set1.write(out);
		
		BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
		
	}
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(Constants)

BOOST_AUTO_TEST_CASE(mapping)
{
	using namespace mapping;
	BOOST_CHECK_EQUAL(chromosome["chr1"], 1);
	BOOST_CHECK_EQUAL(chromosome["chr9"], 9);
	BOOST_CHECK_EQUAL(chromosome["X"], 23);
	BOOST_CHECK_EQUAL(chromosome["chrY"], 24);
	BOOST_CHECK_EQUAL(chromosome["30"], 0);
}

BOOST_AUTO_TEST_SUITE_END()
