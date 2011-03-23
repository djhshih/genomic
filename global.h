#ifndef genomic_global_h
#define genomic_global_h

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdio>

using namespace std;

typedef unsigned long position;

// Number of human autosomes
extern size_t nAutosomes;
// Number of human chromosomes: 22 autosomes + 2 sex chromosomes
extern size_t nChromosomes;


namespace data
{
	enum Type {
		generic, raw, segmented
	};
};

namespace mapping
{
	class ChromosomesMap
	{
	private:
		map<string, size_t> index;
	public:
		ChromosomesMap() {
			// Map numeric values of chromosomes and alternatives prefixed with chr
			char s[6];
			for (size_t i = 1; i <= nChromosomes; ++i) {
				sprintf(s, "%d", i); index[s] = i;
				sprintf(s, "chr%d", i);
				index[s] = i;
			}
			// Map alternative names for chromosomes
			index["X"] = 23;
			index["chrX"] = 23;
			index["Y"] = 24;
			index["chrY"] = 24;
		}	
		size_t operator[] (const string& chr) {
			return index[chr];
		}
		void print() {
			map<string, size_t>::const_iterator end = index.end();
			for (map<string, size_t>::iterator it = index.begin(); it != end; ++it) {
				cout << it->first << " -> " << it->second << endl;
			}
		}
	};
	static ChromosomesMap chromosome;
	
	class ExtensionMap
	{
	private:
		map<string, data::Type> type;
	public:
		ExtensionMap() {
			type["cn"] = data::raw;
			type["seg"] = data::segmented;
		}
		data::Type operator[] (const string& ext) {
			return type[ext];
		}
	};
	static ExtensionMap extension;
};


class Marker
{
public:
	string name;
	size_t chromosome;
	position pos;
	Marker() {}
	Marker(string markerName, size_t markerChromosome, position markerPosition)
	: name(markerName), chromosome(markerChromosome), pos(markerPosition) {}
};

/*
class MarkerSet
{
	friend class RawSampleSet;
public:
	typedef vector<Marker> Chromosome;
	typedef vector<Chromosome> Chromosomes;
	string platform;
private:
	Chromosomes markers;
};
*/

#endif