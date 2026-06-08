#ifndef cna_Marker_h
#define cna_Marker_h

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
#include "Properties.hpp"

namespace marker
{

	class Marker
	{
	public:
		std::string name;
		chromid chromosome;
		position pos;
		bool flag;
		
		Marker() : flag(false) {}
		
		Marker(std::string markerName, chromid markerChromosome, position markerPosition)
		: name(markerName), chromosome(markerChromosome), pos(markerPosition), flag(false) {}
		
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
		typedef std::vector<Marker*> ChromosomeMarkers;
		typedef std::vector<ChromosomeMarkers> GenomeMarkers;
		
		std::string platform;
		bool namedMarkers;
		
		Set(const std::string& markerSetPlatform)
		: platform(markerSetPlatform), namedMarkers(true),
		  set(nChromosomes), unsortedChromIndex(0), refCount(1) {
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
				throw std::logic_error("Unsorted chromosome has yet been initialized.");
			}
		}
		size_t size() const {
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
		
		void setIO(const IOProperties& io) {
			this->io = io;
		}
		
		void read(const std::string& fileName, const std::string& platform, bool doSort=true, bool named=false);
		void read(std::ifstream& file, const std::string& platform, bool doSort=true, bool named=false);
		
		void write(const std::string& fileName);
		void write(std::ofstream& file);
		
		void distribute();
		
		void sort();
		
		bool empty();
		
		// Remove flagged markers from chromosome vectors and delete the owned Marker objects.
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
		
		// Owned marker storage, partitioned by chromosome.
		GenomeMarkers set;
		size_t unsortedChromIndex;
		size_t refCount;
		
		IOProperties io;
		
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
		typedef std::map<std::string, Set*>::iterator iterator;
		std::map<std::string, Set*> sets;
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
		
		void newSetName(std::string& markerSetPlatform);
		
		Set* create(const std::string& markerSetPlatform);
		
		Set* operator[](const std::string& markerSetPlatform) {
			Set* set = sets[markerSetPlatform];
			set->ref();
			return set;
		}
		
		void ref(Set* set) {
			if (set != NULL) {
				ref(set->platform);
			}
		}
		
		void ref(const std::string& markerSetPlatform) {
			iterator it = sets.find(markerSetPlatform);
			if (it != sets.end()) {
				it->second->ref();
			}
		}
		
		void unref(Set* set) {
			if (set != NULL) {
				unref(set->platform);
			}
		}
		
		void unref(const std::string& markerSetPlatform) {
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
