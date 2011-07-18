#include "Sample.h"

const char SampleSet::delim = '\t';

void GenericSampleSet::_read(fstream& file)
{
	string ext = fileName.substr(fileName.find_last_of('.')+1);
	switch (mapping::extension[ext]) {
		case data::raw:
			rep = new RawSampleSet(markers);
			break;
		case data::segmented:
			rep = new SegmentedSampleSet(markers);
			break;
		default:
			throw runtime_error("Cannot determine file type from file name extension");
	}
	rep->_read(file);
}

void GenericSampleSet::_write(fstream& file)
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
	rep->_write(file);
}

RawSampleSet::RawSampleSet(SegmentedSampleSet& set)
{
	//TODO
	// N.B. throw error if set.markers does not exist
	//   Marker information in set is required for creating RawSampleSet
}

void RawSampleSet::_read(fstream& file)
{
	// assume M x (3+N) data matrix with M makers and N samples
	// columns: marker, chromosome, position, samples...
	
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
					create(sampleName);
				}
			} else {
				istringstream stream(line);
				position pos;
				stream >> markerName >> chromName >> pos;
				size_t chr = mapping::chromosome[chromName];
				// ignore unknown chromosome: continue to next line
				if (chr == 0) continue;
				// create marker
				marker::Marker marker(markerName, chr, pos);
				markers->at(chr-1).push_back(marker);
				
				float value;
				size_t i = -1;
				while (!stream.eof()) {
					stream >> value;
					// create point at specified chromosome
					samples[++i]->addToChromosome(chromName, value);
				}
			}
		} else {
			// discard line
		}
	}
}

void RawSampleSet::_write(fstream& file)
{
	file << "marker" << delim << "chromosome" << delim << "position";
	
	// print sample names
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		file << delim << (**it).name;
	}
	file << endl;
	
	// iterate through each chromosome in the vector of vector $markers
	for (size_t chr = 0; chr < markers->size(); ++chr) {
		for (unsigned long markerIndex = 0; markerIndex < markers->at(chr).size(); ++markerIndex) {
			
			// print marker information
			file << markers->at(chr)[markerIndex].name << delim << markers->at(chr)[markerIndex].chromosome << delim << markers->at(chr)[markerIndex].pos;
			
			// iterate through samples to print values, selected the specified chromosome and marker
			SamplesIterator it, end = samples.end();
			for (it = samples.begin(); it != end; ++it) {
				file << delim << (**it)[chr]->at(markerIndex);
			}
			file << endl;
			
		}
	}
}

void RawSampleSet::sort()
{
	// Construct order vector for obtaining a sorted index of markers
	vector< pair<position, size_t> > order;
	
	for (size_t chri = 0; chri < markers->size(); ++chri) {	
		
		// Construct order vector for obtaining a sorted index of markers
		// Additionally, replicate the markers on the curent chromosome;
		//               replicate the chromosome for all samples
		marker::Set::ChromosomeMarkers chromosomeMarkers;
		vector<RawChromosome> samplesChromosomeCopy;
		vector< pair<position, size_t> > order;
		for (size_t j = 0; j < markers->at(chri).size(); ++j) {
			order.push_back( make_pair(markers->at(chri)[j].pos, j) );
			
			chromosomeMarkers.push_back(markers->at(chri)[j]);
			
			for (size_t s = 0; s < samples.size(); ++s) {
				samplesChromosomeCopy.push_back(*(samples[s]->chromosome(chri)));
			}
		}
		
		// Sort on the order vector instead of the original vector<Markers*>,
		//   in order to obtain the sorted index
		std::sort(order.begin(), order.end(), &compare::pair<position, size_t>);
		
		// Now, order is sorted by genomic position, and order[i].second
		//   contains each sorted index
		
		// Sort the markers on current chromosome, and each sample
		for (size_t j = 0; j < markers->at(chri).size(); ++j) {
			size_t index = order[j].second;
			// Set the marker to the corresponding sorted marker
			markers->at(chri)[j] = chromosomeMarkers[index];
			
			// Iterate through samples, set the corresponding marker in the current chromosome
			for (size_t s = 0; s < samples.size(); ++s) {
				samples[s]->chromosome(chri)->at(j) = samplesChromosomeCopy[s][index];
			}
		}
		
	}
}

SegmentedSampleSet::SegmentedSampleSet(RawSampleSet& raw)
{
	clear();
	// use iterators to avoid assuming RawSampleSet stores data in vectors
	// however, need to assume that markers are stored in vectors, for looking up marker information
	
	typedef RawSampleSet::Samples::iterator RawSamplesIterator;
	typedef RawSampleSet::ChromosomesIterator RawChromosomesIterator;
	typedef RawSampleSet::DataIterator RawDataIterator;
	
	// iterate through samples in $raw
	RawSamplesIterator it, end = raw.samples.end();
	for (it = raw.samples.begin(); it != end; ++it) {
		// create sample
		SegmentedSample* sam = create((*it)->name);
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
				prevValue = *markerIt;
				// start at second marker
				++markerIt;
				markerIndex = 1;
				while (markerIt != markerEnd) {
					if (*markerIt != prevValue) {
						// segment ended: store segment from $startMarkerIndex to $markerIndex-1
						Segment seg(
							raw.markers->at(chr)[startMarkerIndex].pos,
							raw.markers->at(chr)[markerIndex-1].pos,
							(markerIndex-1) - startMarkerIndex + 1,
							prevValue
						);
						sam->chromosome(chr)->push_back(seg);
						// start new segment
						startMarkerIndex = markerIndex;
						prevValue = *markerIt;
					}
					++markerIndex;
					++markerIt;
				}
				// store last segment
				// handling is same whether last segment is the last marker alone or
				// 	laste segment ends on the last marker
				Segment seg(
					raw.markers->at(chr)[startMarkerIndex].pos,
					raw.markers->at(chr)[markerIndex-1].pos,
					(markerIndex-1) - startMarkerIndex + 1,
					prevValue
				);
				sam->chromosome(chr)->push_back(seg);
			}
			++chr;
		}  // for chrIt
	} // for it
}

void SegmentedSampleSet::_read(fstream& file)
{
	// assume M x 6 data matrix
	// columns: sample, chr, start, end, markers, value
	
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
			Segment seg;
			file >> seg.start >> seg.end >> seg.nelem >> seg.value;
			create(sampleName)->addToChromosome(chromName, seg);
			//trace("%s %s %d %d %d %f\n", sampleName.c_str(), chromName.c_str(), seg->start, seg->end, seg->nelem, seg->value);
		} else {
			// discard line
			getline(file, line);
		}
	}
}

void SegmentedSampleSet::_write(fstream& file)
{
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
				file << (*it)->name << delim << chr << delim << segIt->start << delim << segIt->end << delim << segIt->nelem << delim << segIt->value << endl;
			}
			++chr;
		}
	}
}

void SegmentedSampleSet::sort()
{
	// Sort samples by name
	std::sort(samples.begin(), samples.end(), &SegmentedSample::pcompare);
	// Iterate through samples and chromosomes therefore, sort segments
	SamplesIterator it, end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		ChromosomesIterator chrIt, chrEnd = (*it)->end();
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			std::sort(chrIt->begin(), chrIt->end(), &Segment::compare);
		}
	}
}

void SegmentedSampleSet::filter(SegmentedSampleSet& ref)
{
	// iteratrate through samples
	SamplesIterator it;
	const SamplesIterator end = samples.end();
	for (it = samples.begin(); it != end; ++it) {
		// iterate through chromosomes
		ChromosomesIterator chrIt;
		const ChromosomesIterator chrEnd = (*it)->end();
		size_t chri = 0;
		for (chrIt = (*it)->begin(); chrIt != chrEnd; ++chrIt) {
			// iterate through segments on a chromosome
			DataIterator segIt;
			const DataIterator segEnd = chrIt->end();
			for (segIt = chrIt->begin(); segIt != segEnd; ++segIt) {
				// Compare against segments on reference
				//TODO
			}
			++chri;
		}
	}
}
