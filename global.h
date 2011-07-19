#ifndef genomic_global_h
#define genomic_global_h

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <cmath>
#include <string>

using namespace std;

typedef unsigned long position;

// Number of human autosomes
extern size_t nAutosomes;
// Number of human chromosomes: 22 autosomes + 2 sex chromosomes
extern size_t nChromosomes;


namespace data
{
	enum Type {
		generic, raw, segmented, raw_ascn, segmented_ascn, dchip, cnag, picnic
	};
};

namespace mapping
{
	/* singleton */
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
	
	/* singleton */
	class ExtensionMap
	{
	private:
		map<string, data::Type> type;
	public:
		ExtensionMap() {
			type["cn"] = data::raw;
			type["seg"] = data::segmented;
			type["ascn"] = data::raw_ascn;
			type["asseg"] = data::segmented_ascn;
			type["dchip"] = data::dchip;
			type["cnag"] = data::cnag;
			type["picnic"] = data::picnic;
		}
		data::Type operator[] (const string& ext) {
			return type[ext];
		}
	};
	static ExtensionMap extension;
};

namespace compare
{
	template <typename T1, typename T2>
	bool pair(pair<T1, T2> a, pair<T1, T2> b) {
		return a.first < b.first;
	}
}

namespace marker
{

	class Marker
	{
	public:
		string name;
		size_t chromosome;
		position pos;
		Marker() {}
		Marker(string markerName, size_t markerChromosome, position markerPosition)
		: name(markerName), chromosome(markerChromosome), pos(markerPosition) {}
		static bool compare(const Marker& a, const Marker& b) {
			return a.pos < b.pos;
		}
		static bool pcompare(Marker* a, Marker* b) {
			return a->pos < b->pos;
		}
	};

	class Set
	{
		friend class Manager;
	public:
		typedef vector<Marker> ChromosomeMarkers;
		typedef vector<ChromosomeMarkers> GenomeMarkers;
		typedef typename GenomeMarkers::iterator ChromosomesIterator;
		string platform;
		Set(const string& markerSetPlatform) : refCount(1), set(nChromosomes), platform(markerSetPlatform) {
		}
		ChromosomeMarkers& at(size_t i) {
			return set[i];
		}
		ChromosomeMarkers& operator[](size_t i) {
			return set[i];
		}
		const size_t size() const {
			return set.size();
		}
		ChromosomesIterator begin() {
			return set.begin();
		}
		ChromosomesIterator end() {
			return set.end();
		}
		void read(const string& fileName, const string& platform) {
			ifstream file(fileName.c_str(), ios::in);
			if (!file.is_open()) throw runtime_error("Failed to open input file");
			read(file, platform);
			file.close();
		}
		void read(ifstream& file, const string& platform) {
			clear();
			// assume M x 3 data matrix with M makers
			// columns: marker, chromosome, position
			
			string line;
			size_t lineCount = 0, nSkippedLines = 0, headerLine = 1;
			string markerName, chromName, sampleName, discard;
			while (true) {
				getline(file, line);
				
				if (file.eof()) break;
				if (++lineCount > nSkippedLines) {
					if (lineCount == headerLine) {
						// discard
					} else {
						istringstream stream(line);
						position pos;
						stream >> markerName >> chromName >> pos;
						size_t chr = mapping::chromosome[chromName];
						// ignore unknown chromosome: continue to next line
						if (chr == 0) continue;
						// create marker
						marker::Marker marker(markerName, chr, pos);
						set[chr-1].push_back(marker);
					}
				} else {
					// discard line
				}
			}
			
			sort();
		}
		
		void sort() {
			ChromosomesIterator it;
			const ChromosomesIterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				std::sort(it->begin(), it->end(), &Marker::compare);
			}
		}
		
		bool empty() {
			bool isEmpty = true;
			GenomeMarkers::iterator it;
			const GenomeMarkers::iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				isEmpty &= it->empty();
			}
			return isEmpty;
		}
		
	private:
		GenomeMarkers set;
		size_t refCount;
		void ref() {
			++refCount;
		}
		bool unref() {
			if (--refCount < 1) {
				clear();
				return true;
			}
			return false;
		}
		void clear() {
			GenomeMarkers::iterator it;
			const GenomeMarkers::iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				it->clear();
			}
		}
	};

	/* singleton */
	// Manage marker sets
	class Manager
	{
	private:
		typedef map<string, Set*>::iterator iterator;
		map<string, Set*> sets;
	public:
		Manager() {}
		~Manager() {
			iterator it;
			const iterator end = sets.end();
			for (it = sets.begin(); it != end; ++it) {
				delete it->second;
			}
		}
		Set* create(const string& markerSetPlatform) {
			iterator it = sets.find(markerSetPlatform);
			Set* set;
			if (it == sets.end()) {
				set = new Set(markerSetPlatform);
				sets[markerSetPlatform] = set;
			} else {
				set = it->second;
				set->ref();
			}
			return set;
		}
		Set* operator[](const string& markerSetPlatform) {
			Set* set = sets[markerSetPlatform];
			set->ref();
			return set;
		}
		void unref(Set* set) {
			if (set != NULL) unref(set->platform);
		}
		void unref(const string& markerSetPlatform) {
			iterator it = sets.find(markerSetPlatform);
			if (it != sets.end() && it->second->unref()) {
				delete (it->second);
				sets.erase(it);
			}
		}
	};
	static Manager manager;

}

namespace name
{
	string common(const string& a, const string& b);
}

#endif