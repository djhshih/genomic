#include "Marker.hpp"

namespace marker {

	void Set::read(const string& fileName, const string& platform, bool doSort) {
		ifstream file(fileName.c_str(), ios::in);
		if (!file.is_open()) throw runtime_error("Failed to open input file");
		read(file, platform, doSort);
		file.close();
	}
	
	void Set::read(ifstream& file, const string& platform, bool doSort) {
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
	
	void Set::sort() {
		GenomeMarkers::iterator it;
		GenomeMarkers::const_iterator end = set.end();
		for (it = set.begin(); it != end; ++it) {
			std::sort(it->begin(), it->end(), &Marker::pcompare);
		}
	}
	
	bool Set::empty() {
		bool isEmpty = true;
		GenomeMarkers::iterator it;
		GenomeMarkers::const_iterator end = set.end();
		for (it = set.begin(); it != end; ++it) {
			isEmpty &= it->empty();
		}
		return isEmpty;
	}
	
	void Set::distribute() {
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
	
	// TODO make more efficient
	void Set::clean() {
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
	
	// construct hash containing all marker names found in markerNames
	void Set::hashMarkers(const std::vector<std::string>& markerNames, uset& hash) {
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
	void Set::hashMarkers(const Set& ref, uset& hash) {
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
	
	void Set::filter(const uset& refMarkersHash) {
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
	
	void Set::clear() {
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

	const char Manager::alphanum[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	// length of alphanum, subtracted by 1 (to account for null character)
	const unsigned Manager::nalphanum = 62;
	
	void Manager::newSetName(string& markerSetPlatform) {
		const unsigned n = 7;
		char randstr[n];
		do {
			randAlphaNum(randstr, n);
			markerSetPlatform = randstr;
		} while (sets.find(markerSetPlatform) != sets.end());
	}

	Set* Manager::create(const string& markerSetPlatform) {
		if (markerSetPlatform.empty()) {
			throw runtime_error("Market set platform name cannot be empty");
		}
		Set* set;
		iterator it = sets.find(markerSetPlatform);
		if (it == sets.end()) {
			set = new Set(markerSetPlatform);
			sets[markerSetPlatform] = set;
		} else {
			set = it->second;
			set->ref();
		}
		return set;
	}
	
}
