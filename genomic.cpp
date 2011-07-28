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
			fileNames.push_back(name::filepath(samplesFileName) + sample);
			cout << sample << endl;
		}
	}
	
	PicnicSampleSet pset;
	pset.read(fileNames, markersFileName, true);
	pset.write("output.picnic.ascn");
	
	SegmentedSampleSet<PicnicSampleSet::Value> sset(pset);
	sset.write("output.asseg");

	return 0;
}
