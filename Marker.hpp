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

#include <boost/unordered_set.hpp>

#include "global.hpp"

namespace marker
{

	class Marker
	{
	public:
		string name;
		chromid chromosome;
		position pos;
		bool flag;
		
		Marker() : flag(false) {}
		
		Marker(string markerName, chromid markerChromosome, position markerPosition)
		: flag(false), name(markerName), chromosome(markerChromosome), pos(markerPosition) {}
		
		static bool compare(const Marker& a, const Marker& b) {
			return a.pos < b.pos;
		}
		
		static bool pcompare(Marker* a, Marker* b) {
			return a->pos < b->pos;
		}
	};

	//TODO store Marker as pointers
	class Set
	{
		friend class Manager;
		typedef boost::unordered_set<std::string> uset;
	public:
		typedef vector<Marker*> ChromosomeMarkers;
		typedef vector<ChromosomeMarkers> GenomeMarkers;
		
		string platform;
		
		Set(const string& markerSetPlatform)
		: refCount(1), set(nChromosomes), unsortedChromIndex(0),
			platform(markerSetPlatform) {
		}
		
		~Set() {
			clear();
		}
		
		ChromosomeMarkers& at(size_t i) {
			return set[i];
		}
		ChromosomeMarkers& operator[](size_t i) {
			return set[i];
		}
		ChromosomeMarkers& unsortedChromosome() {
			if (unsortedChromIndex > 0) {
				return set[unsortedChromIndex];
			} else {
				throw logic_error("Unsorted chromosome has yet been initialized.");
			}
		}
		const size_t size() const {
			return set.size();
		}
		GenomeMarkers::const_iterator begin() {
			return set.begin();
		}
		GenomeMarkers::const_iterator end() {
			return set.end();
		}
		
		void addToChromosome(size_t chromIndex, Marker* marker) {
			set[chromIndex].push_back(marker);
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
						marker::Marker* marker = new Marker(markerName, chr, pos);
						
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
				ChromosomeMarkers::const_iterator it, end = set[unsortedChromIndex].end();
				for (it = set[unsortedChromIndex].begin(); it != end; ++it) {
					set[(*it)->chromosome-1].push_back(*it);
				}
				// Remove extra chromosome
				set[unsortedChromIndex].clear();
				set.resize(set.size()-1);
				unsortedChromIndex = 0;
			}
		}
		
		void sort() {
			GenomeMarkers::iterator it;
			GenomeMarkers::const_iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				std::sort(it->begin(), it->end(), &Marker::pcompare);
			}
		}
		
		bool empty() {
			bool isEmpty = true;
			GenomeMarkers::iterator it;
			GenomeMarkers::const_iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				isEmpty &= it->empty();
			}
			return isEmpty;
		}
		
		// TODO make more efficient
		void clean() {
			GenomeMarkers::iterator it;
			GenomeMarkers::const_iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				ChromosomeMarkers::iterator chromIt = it->begin();
				while (chromIt != it->end()) {
					// VERY INEFFICIENT!
					if ((*chromIt)->flag) {
						//trace("Erasing %s...\n", chromIt->name.c_str());
						it->erase(chromIt);
					} else {
						++chromIt;
					}
				}
			}
		}
		
		void filter(const std::vector<std::string>& refMarkers) {
			uset hash;
			hashMarkers(refMarkers, hash);
			filter(hash);
		}
		
		// flag markers found in ref
		void filter(const Set& ref) {
			uset hash;
			hashMarkers(ref, hash);
			filter(hash);
		}
		
	private:
		
		GenomeMarkers set;
		size_t unsortedChromIndex;
		size_t refCount;
		
		char delim;
		size_t nSkippedLines;
		size_t headerLine;
		
		// construct hash containing all marker names found in markerNames
		void hashMarkers(const std::vector<std::string>& markerNames, uset& hash) {
			if (markerNames.size() == 0) return;
			// initialize hash with at least enough buckets to hold all markers,
			//  multipled by factor to for expected load factor of 0.7
			hash.rehash(markerNames.size() * 1.5);
			
			// populate hash
			std::vector<std::string>::const_iterator it, end = markerNames.end();
			for (it = markerNames.begin(); it != end; ++it) {
				hash.insert(*it);
			}
		}
		
		// construct hash containing all marker names found in ref
		void hashMarkers(const Set& ref, uset& hash) {
			// determine total number of reference markers
			GenomeMarkers::const_iterator end = ref.set.end();
			size_t numRefMarkers = 0;
			for (GenomeMarkers::const_iterator it = ref.set.begin(); it != end; ++it) {
				numRefMarkers += it->size();
			}
			
			// initialize hash with at least enough buckets to hold all markers,
			//  multipled by factor to for expected load factor of 0.7
			hash.rehash(numRefMarkers * 1.5);
			
			// populate hash
			for (GenomeMarkers::const_iterator it = ref.set.begin(); it != end; ++it) {
				ChromosomeMarkers::const_iterator markerIt, markerEnd = it->end();
				for (markerIt = it->begin(); markerIt != markerEnd; ++markerIt) {
					hash.insert((*markerIt)->name);
				}
			}
		}
		
		void filter(const uset& refMarkersHash) {
			// flag markers in this->set that are found in refMarkers
			GenomeMarkers::const_iterator end = set.end();
			size_t filterCount = 0;
			for (GenomeMarkers::iterator it = set.begin(); it != end; ++it) {
				ChromosomeMarkers::iterator markerIt;
				ChromosomeMarkers::const_iterator markerEnd = it->end();
				uset::const_iterator refMarkerEnd = refMarkersHash.end();
				for (markerIt = it->begin(); markerIt != markerEnd; ++markerIt) {
					if (refMarkersHash.find((*markerIt)->name) != refMarkerEnd) {
						(*markerIt)->flag = true;
						++filterCount;
					}
				}
			}
			
			trace("Number of markers filtered: %d\n", filterCount);
		}
		
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
			GenomeMarkers::const_iterator end = set.end();
			for (it = set.begin(); it != end; ++it) {
				ChromosomeMarkers::iterator markerIt;
				ChromosomeMarkers::const_iterator markerEnd = it->end();
				for (markerIt = it->begin(); markerIt != markerEnd; ++markerIt) {
					delete (*markerIt);
				}
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
