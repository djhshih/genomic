#include "genomic.h"
#include "Sample.h"

#include <cstdlib>
#include <fstream>

using namespace std;


int main(int argc, char **argv)
{
	
	string markersFileName = argv[1];
	string samplesFileName = argv[2];

	ifstream samplesFile(samplesFileName.c_str());
	if (!samplesFile.is_open()) {
		throw runtime_error("Failed to open samples file");
	}
	
	vector<string> fileNames;
	while (!samplesFile.eof()) {
		string sample;
		samplesFile >> sample;
		if (sample != "") {
			fileNames.push_back(sample);
			cout << sample << endl;
		}
	}
	
	//vector<string> fileNames(2);
	//fileNames[0] = "tests/picnic1a.in";
	//fileNames[1] = "tests/picnic1b.in";
	//string markersFileName = "tests/picnic.snp6.tsv";
	//string markersFileName = "tests/picnic.snp6.csv";
	
	PicnicSampleSet pset;
	pset.read(fileNames, markersFileName);
	//pset.write("tests/picnic1.out");
	pset.write("output.picnic.ascn");
	
	SegmentedSampleSet<PicnicSampleSet::Value> sset(pset);
	//sset.write("tests/picnic1.seg");
	sset.write("output.asseg");

	return 0;
}
