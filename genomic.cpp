#include "genomic.h"

#include "Tree.h"
#include <cstdlib>

using namespace std;

class Datum {
	int* rep;
	int ndim;
public:
	Datum(size_t _ndim) : ndim(_ndim), rep(new int[ndim]) {}
	Datum(const Datum& d) : ndim(d.ndim), rep(new int[d.ndim]) {
		for (size_t i = 0; i < d.ndim; ++i) {
			rep[i] = d.rep[i];
		}
	}
	~Datum() {
		delete rep;
	}
	Datum& operator=(const Datum& d) {
		Datum tmp(d);
		swap(this->rep, tmp.rep);
		return *this;
	}
	int& operator[](size_t i) {
		return rep[i];
	}
	int* keys() const {
		return rep;
	}
};

int main(int argc, char **argv)
{
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
	tree::KDTree<Datum> tree2(ndim);
	for (int i = 0; i < 20; ++i) {
		Datum dat(ndim);
		for (int d = 0; d < ndim; ++d) {
			dat[d] = (rand() % 100) + 1;
		}
		/*
		// print tuple
		cout << "/";
		for (int i = 0; i < ndim; ++i) {
			cout << dat.keys()[i] << "/";
		}
		cout << endl;
		*/
		tree2.insert(dat);
	}
	tree2.print();
	

	return 0;
}
