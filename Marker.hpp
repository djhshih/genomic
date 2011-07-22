
#ifndef genomic_Marker_h
#define genomic_Marker_h

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <algorithm>

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
		typedef typename ChromosomeMarkers::iterator MarkersIterator;
		string platform;
		Set(const string& markerSetPlatform)
		: refCount(1), set(nChromosomes), unsortedChromIndex(0),
			platform(markerSetPlatform) {
		}
		ChromosomeMarkers& at(size_t i) {
			return set[i];
		}
		ChromosomeMarkers& operator[](size_t i) {
			return set[i];
		}
		ChromosomeMarkers* unsortedChromosome() {
			if (unsortedChromIndex > 0) {
				return &set[unsortedChromIndex];
			}
			return NULL;
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
		void setIO(char _delim, size_t _headerLine, size_t _nSkippedLines) {
			delim = _delim;
			headerLine = _headerLine;
			nSkippedLines = _nSkippedLines;
		}
		void read(const string& fileName, const string& platform, bool doSort=true) {
			ifstream file(fileName.c_str(), ios::in);
			if (!file.is_open()) throw runtime_error("Failed to open input file");
			read(file, platform, doSort);
			file.close();
		}
		void read(ifstream& file, const string& platform, bool doSort=true) {
			clear();
			// assume M x 3 data matrix with M makers
			// columns: marker, chromosome, position
			
			if (!doSort) {
				unsortedChromIndex = set.size();
				set.resize(set.size()+1);
			}
			
			size_t lineCount = 0;
			string line;
			position pos;
			string markerName, chromName, s;
			while (true) {
				getline(file, line);
				
				if (file.eof()) break;
				if (++lineCount > nSkippedLines && lineCount != headerLine) {
						istringstream stream(line);
						//stream >> markerName >> chromName >> pos;
						getline(stream, markerName, delim);
						getline(stream, chromName, delim);
						getline(stream, s, delim);
						pos = atol(s.c_str());
						
						size_t chr = mapping::chromosome[chromName];
						// ignore unknown chromosome: continue to next line
						if (chr == 0) continue;
						// create marker
						marker::Marker marker(markerName, chr, pos);
						
						if (doSort) {
							// add marker onto appropriate chromosome
							set[chr-1].push_back(marker);
						} else {
							// add all markers to extra chromosome
							set[unsortedChromIndex].push_back(marker);
						}
				} else {
					// discard line
				}
			}
			
			if (doSort) sort();
		}
		
		void distribute() {
			if (unsortedChromIndex > 0) {
				// Move markers from the first chromosome onto appropriate chromosomes
				MarkersIterator it;
				const MarkersIterator end = set[unsortedChromIndex].end();
				for (it = set[unsortedChromIndex].begin(); it != end; ++it) {
					set[it->chromosome-1].push_back(*it);
				}
				// Remove extra chromosome
				set[unsortedChromIndex].clear();
				set.resize(set.size()-1);
				unsortedChromIndex = 0;
			}
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
		size_t unsortedChromIndex;
		size_t refCount;
		
		char delim;
		size_t nSkippedLines;
		size_t headerLine;
		
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
			unsortedChromIndex = 0;
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

#endif
