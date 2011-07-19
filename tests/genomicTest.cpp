#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Unit Tests for genomic"
#include <boost/test/unit_test.hpp>

#include "../global.h"
#include "../Sample.h"
#include "FilesDiff.hpp"

#include <fstream>
#include <queue>
#include <stdexcept>

using namespace std;

template <typename CopyNumberValueType> struct IOFixture;

template <typename CopyNumberValueType>
class IOTests {
	friend struct IOFixture<CopyNumberValueType>;
public:
	enum TestType { Basic = 0, Conversion };
private:
	vector<string> configFilenames;
	IOTests() {
		configFilenames.push_back("iotest_basic.cfg");
	}
};

template <>
class IOTests<AlleleSpecificCopyNumberValue> {
	friend struct IOFixture<AlleleSpecificCopyNumberValue>;
private:
	vector<string> configFilenames;
	IOTests() {
		configFilenames.push_back("iotest_basic_as.cfg");
	}
};

template <typename CopyNumberValueType>
struct IOFixture
{
	FilesDiff diff;
	fstream configf;
	IOTests<CopyNumberValueType> tests;
	const char chrComment, chrClass;
	
	RawSampleSet<CopyNumberValueType> rset;
	SegmentedSampleSet<CopyNumberValueType> sset;
	GenericSampleSet<CopyNumberValueType> gset;
	
	vector< queue<SampleSet<CopyNumberValueType>*> > sets;
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
					} else if (s == "SegmentedSampleSet") {
						sets[type].push(&sset);
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

BOOST_FIXTURE_TEST_CASE(InputOutput, IOFixture<CopyNumberValue>)
{
	size_t type;
	
	BOOST_TEST_MESSAGE("Basic IO");
	type = (size_t)IOTests<CopyNumberValue>::Basic;
	while (!sets[type].empty()) {
		SampleSet<CopyNumberValue>* set = sets[type].front();
		sets[type].pop();
		set->read(filenames[type].front());
		filenames[type].pop();
		string out = filenames[type].front();
		filenames[type].pop();
		set->write(out);
		string ans = filenames[type].front();
		filenames[type].pop();
		BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
	}
	if (!filenames[type].empty()) {
		throw runtime_error("filenames queue was not exhausted: SampleSet IO Test configuration file is likely malformed");
	}
}

BOOST_FIXTURE_TEST_CASE(InputOutput_AlleleSpecific, IOFixture<AlleleSpecificCopyNumberValue>)
{
	size_t type;
	
	BOOST_TEST_MESSAGE("Basic IO (Allele-specific)");
	type = (size_t)IOTests<CopyNumberValue>::Basic;
	while (!sets[type].empty()) {
		SampleSet<AlleleSpecificCopyNumberValue>* set = sets[type].front();
		sets[type].pop();
		set->read(filenames[type].front());
		filenames[type].pop();
		string out = filenames[type].front();
		filenames[type].pop();
		set->write(out);
		string ans = filenames[type].front();
		filenames[type].pop();
		BOOST_CHECK_EQUAL(diff.different(out, ans), 0);
	}
	if (!filenames[type].empty()) {
		throw runtime_error("filenames queue was not exhausted: SampleSet IO Test configuration file is likely malformed");
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