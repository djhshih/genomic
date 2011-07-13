#include "Sample.h"

const char SampleSet::delim = '\t';

void GenericSampleSet::read(const string& fileName)
{
	clear();
	string ext = fileName.substr(fileName.find_last_of('.')+1);
	switch (mapping::extension[ext]) {
		case data::raw:
			rep = new RawSampleSet();
			break;
		case data::segmented:
			rep = new SegmentedSampleSet();
			break;
		default:
			throw runtime_error("Cannot determine file type from file name extension");
	}
	rep->read(fileName);
}

void GenericSampleSet::write(const string& fileName)
{
	string ext = fileName.substr(fileName.find_last_of('.')+1);
	
	// cast $rep to appropriate type
	// runtime checking is skipped (i.e. use static_cast instead of dynamic_case)
	//  since exact type can be determined
	SampleSet* tmp = NULL;
	switch (mapping::extension[ext]) {
		case data::raw:
			switch (rep->type()) {
				case data::raw:
					// nothing needs to be done
					break;
				case data::segmented:
					tmp = new RawSampleSet(*static_cast<SegmentedSampleSet*>(rep));
					swap(rep, tmp);
					delete tmp;
					break;
			}
			break;
		case data::segmented:
			switch (rep->type()) {
				case data::raw:
					tmp = new SegmentedSampleSet(*static_cast<RawSampleSet*>(rep));
					swap(rep, tmp);
					delete tmp;
					break;
				case data::segmented:
					// nothing to be done
					break;
			}
			break;
		default:
			throw runtime_error("Cannot determine file type from file name extension");
	}
	rep->write(fileName);
}

RawSampleSet::RawSampleSet(SegmentedSampleSet& set)
{
	//TODO
}

void RawSampleSet::read(const string& fileName)
{
	clear();
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	file.open(fileName.c_str(), ios::in);
	if (!file.is_open()) throw runtime_error("Failed to open input file");
	
	string line;
	
	size_t nSkippedLines = 0;
	size_t headerLine = 1;
	size_t lineCount = 0;
	string markerName, chromName, sampleName, discard;
	while (true) {
		getline(file, line);
		
		if (file.eof()) break;
		if (++lineCount > nSkippedLines) {
			if (lineCount == headerLine) {
				istringstream stream(line);
				// discard the marker information columns (3)
				stream >> discard >> discard >> discard;
				// create samples
				while (!stream.eof()) {
					stream >> sampleName;
					// create sample with $sampleName	
					sample(sampleName);
				}
			} else {
				istringstream stream(line);
				position pos;
				stream >> markerName >> chromName >> pos;
				size_t chr = mapping::chromosome[chromName];
				// ignore unknown chromosome: continue to next line
				if (chr == 0) continue;
				// create marker
				markers[chr-1].push_back(new Marker(markerName, chr, pos));
				
				float value;
				size_t i = -1;
				while (!stream.eof()) {
					stream >> value;
					// create point at specified chromosome
					Value* point = samples[++i]->createAt(chromName);
					if (point == NULL) {
						throw runtime_error("Failed to create point on specified chromosome");
					}	
					*point = value;
				}
			}
		} else {
			// discard line
		}
	}
	
	file.close();
	trace("Read file %s\n", fileName.c_str());
}

void RawSampleSet::write(const string& fileName)
{
	file.open(fileName.c_str(), ios::out);
	if (!file.is_open()) throw runtime_error("Failed to open output file");
	
	file << "marker" << delim << "chromosome" << delim << "position";
	
	// print sample names
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << endl;
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers.size(); ++chr) {
		for (unsigned long markerIndex = 0; markerIndex < markers[chr].size(); ++markerIndex) {
			
			// print marker information
			file << markers[chr][markerIndex]->name << delim << markers[chr][markerIndex]->chromosome << delim << markers[chr][markerIndex]->pos;
			
			// iterate through samples to print values, selected the specified chromosome and marker
			SamplesIterator it, end = samples.end();
			for (it = samples.begin(); it != end; ++it) {
				file << delim << *((**it)[chr]->at(markerIndex));
			}
			file << endl;
			
		}
	}
	
	file.close();
	trace("Wrote file %s\n", fileName.c_str());
}

SegmentedSampleSet::SegmentedSampleSet(RawSampleSet& raw)
{
	clear();
	// use iterators to avoid assuming RawSampleSet stores data in vectors
	// however, need to assume that markers are stored in vectors, for looking up marker information
	
	typedef RawSampleSet::Samples::iterator RawSamplesIterator;
	typedef Sample<RawSampleSet::Value>::Chromosomes::iterator RawChromosomesIterator;
	typedef Sample<RawSampleSet::Value>::Chromosome::iterator RawDataIterator;
	
	// iterate through samples in $raw
	RawSamplesIterator it, end = raw.samples.end();
	for (it = raw.samples.begin(); it != end; ++it) {
		// create sample
		Sample<Segment>* sam = sample((*it)->name);
		// iterate through chromosome in sample
		size_t chr = 0;
		RawChromosomesIterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through markers on chromosome
			// CANNOT assume markers are stored in a vector
			RawDataIterator markerIt = chrIt->begin(), markerEnd = chrIt->end();
			position markerIndex, startMarkerIndex = 0;
			RawSampleSet::Value prevValue;
			if (markerIt != markerEnd) {
				prevValue = **markerIt;
				// start at second marker
				++markerIt;
				markerIndex = 1;
				while (markerIt != markerEnd) {
					if (**markerIt != prevValue) {
						// segment ended: store segment from $startMarkerIndex to $markerIndex-1
						sam->chromosome(chr)->push_back( new Segment(
							raw.markers[chr][startMarkerIndex]->pos,
							raw.markers[chr][markerIndex-1]->pos,
							(markerIndex-1) - startMarkerIndex + 1,
							prevValue) );
						// start new segment
						startMarkerIndex = markerIndex;
						prevValue = **markerIt;
					}
					++markerIndex;
					++markerIt;
				}
				// store last segment
				// handling is same whether last segment is the last marker alone or
				// 	laste segment ends on the last marker
				sam->chromosome(chr)->push_back( new Segment(
					raw.markers[chr][startMarkerIndex]->pos,
					raw.markers[chr][markerIndex-1]->pos,
					(markerIndex-1) - startMarkerIndex + 1,
					prevValue) );
			}
			++chr;
		}  // for chrIt
	} // for it
}

void SegmentedSampleSet::read(const string& fileName)
{
	clear();
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	file.open(fileName.c_str(), ios::in);
	if (!file.is_open()) throw runtime_error("Failed to open input file");
	
	string line;
	size_t nSkippedLines = 1;
	size_t lineCount = 0;
	string sampleName, chromName;
	Segment* seg;
	while (true) {
		if (++lineCount > nSkippedLines) {
			file >> sampleName >> chromName;
			if (file.eof()) break;
			// ignore unknown chromosome: continue to next line
			if (mapping::chromosome[chromName] == 0) continue;
			// create segment at specified chromosome
			seg = sample(sampleName)->createAt(chromName);
			if (seg == NULL) {
				throw runtime_error("Failed to create segment on specified chromosome");
			}	
			file >> seg->start >> seg->end >> seg->nelem >> seg->value;
			//trace("%s %s %d %d %d %f\n", sampleName.c_str(), chromName.c_str(), seg->start, seg->end, seg->nelem, seg->value);
		} else {
			// discard line
			getline(file, line);
		}
	}
	
	file.close();
	trace("Read file %s\n", fileName.c_str());
}

void SegmentedSampleSet::write(const string& fileName)
{
	file.open(fileName.c_str(), ios::out);
	if (!file.is_open()) throw runtime_error("Failed to open output file");
	
	file << "sample" << delim << "chromosome" << delim << "start" << delim << "end" << delim << "count" << delim << "state" << endl;
	
	// iteratrate through samples
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		// can assume that Sample are stored chromosomes in a vector
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		size_t chr = 1;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			// CANNOT assume segments are stored in a vector
			DataIterator segIt, segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				Segment* seg = (*segIt);
				file << (*it)->name << delim << chr << delim << seg->start << delim << seg->end << delim << seg->nelem << delim << seg->value << endl;
			}
			++chr;
		}
	}
	file.close();
	trace("Wrote file %s\n", fileName.c_str());
}