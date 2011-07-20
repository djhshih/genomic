#include "genomic.h"

#include "Tree.h"
#include "Sample.h"

#include <cstdlib>

using namespace std;


int main(int argc, char **argv)
{
	/*
	srand((unsigned int)time(NULL));
	cout << "BSTree" << endl;
	tree::BSTree<int> tree;
	tree.clear();
	for (int i = 0; i < 20; ++i) {
		tree.insert((rand() % 100) + 1);
	}
	tree.print();
	cout << endl;
	
	cout << "KDTree" << endl;
	const int ndim = 3;
	tree::KDTree<tree::Container<int> > tree2(ndim);
	for (int i = 0; i < 20; ++i) {
		tree::Container<int> dat(ndim);
		for (int d = 0; d < ndim; ++d) {
			dat[d] = (rand() % 100) + 1;
		}
		tree2.insert(dat);
	}
	tree2.print();
	*/
	
	vector<string> fileNames(2);
	fileNames[0] = "tests/picnic1a.in";
	fileNames[1] = "tests/picnic1b.in";
	//string markersFileName = "tests/picnic.snp6.tsv";
	string markersFileName = "tests/picnic.snp6.csv";
	
	PicnicSampleSet pset;
	pset.read(fileNames, markersFileName);
	pset.write("tests/picnic1.out");
	
	SegmentedSampleSet<AlleleSpecificCopyNumberValue> sset(pset);
	sset.write("tests/picnic1.seg");

	return 0;
}
