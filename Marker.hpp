#ifndef genomic_Marker_h
#define genomic_Marker_h

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include <boost/unordered_set.hpp>

#include "typedefs.h"
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
		
		ChromosomeMarkers& at(chromid i) {
			return set[i];
		}
		ChromosomeMarkers& operator[](chromid i) {
			return set[i];
		}
		void copyChromosome(chromid chromIndex, ChromosomeMarkers& dest) {
			dest.clear();
			dest.reserve(set[chromIndex].size());
			ChromosomeMarkers::const_iterator it, end = set[chromIndex].end();
			for (it = set[chromIndex].begin(); it != end; ++it) {
				dest.push_back(*it);
			}
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
		
		void addToChromosome(chromid chromIndex, Marker* marker) {
			set[chromIndex].push_back(marker);
		}
		
		void setIO(char _delim, size_t _headerLine, size_t _nSkippedLines) {
			delim = _delim;
			headerLine = _headerLine;
			nSkippedLines = _nSkippedLines;
		}
		
		void read(const string& fileName, const string& platform, bool doSort=true);
		void read(ifstream& file, const string& platform, bool doSort=true);
		
		void distribute();
		
		void sort();
		
		bool empty();
		
		void clean();
		
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
		void hashMarkers(const std::vector<std::string>& markerNames, uset& hash);
		
		// construct hash containing all marker names found in ref
		void hashMarkers(const Set& ref, uset& hash);
		
		void filter(const uset& refMarkersHash);
		
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
		
		void clear();
		
	};

	/* singleton */
	// Manage marker sets
	class Manager
	{
	private:
		typedef map<string, Set*>::iterator iterator;
		map<string, Set*> sets;
		static const char alphanum[];
		static const unsigned nalphanum;
	public:
		Manager() {}
		~Manager() {
			iterator it;
			const iterator end = sets.end();
			for (it = sets.begin(); it != end; ++it) {
				delete it->second;
			}
		}
		
		void newSetName(string& markerSetPlatform);
		
		Set* create(const string& markerSetPlatform);
		
		Set* operator[](const string& markerSetPlatform) {
			Set* set = sets[markerSetPlatform];
			set->ref();
			return set;
		}
		
		void ref(Set* set) {
			if (set != NULL) ref(set->platform);
		}
		
		void ref(const string& markerSetPlatform) {
			iterator it = sets.find(markerSetPlatform);
			if (it != sets.end()) {
				it->second->ref();
			}
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
		
	private:
		
		void randAlphaNum(char *s, const unsigned n) {
			for (size_t i = 0; i < n-1; ++i) {
				s[i] = alphanum[std::rand() % (nalphanum)];
			}
			s[n] = '\0';
		}
		
	};
	static Manager manager;

}

#endif
