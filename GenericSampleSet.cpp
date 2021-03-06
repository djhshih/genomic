#include "GenericSampleSet.hpp"

void GenericSampleSet::_read(fstream& file)
{
	const string& fileName = Base::fileName;
	marker::Set* markers = Base::markers;
	
	string ext = name::fileext(fileName);
	switch (mapping::extension[ext]) {
		case data::raw:
			rep = new RawSampleSet<cnvalue>(markers);
			break;
		case data::raw_ascn:
			rep = new RawSampleSet<alleles_cn>(markers);
			break;
		case data::segmented:
			rep = new SegmentedSampleSet<cnvalue>(markers);
			break;
		case data::segmented_ascn:
			rep = new SegmentedSampleSet<alleles_cn>(markers);
		default:
			throw invalid_conversion("Cannot determine file type from file name extension.");
	}
	rep->_read(file);
}

void GenericSampleSet::_write(fstream& file)
{
	const string& fileName = Base::fileName;
	
	string ext = name::fileext(fileName);
	
	// cast $rep to appropriate type
	// runtime checking is skipped (i.e. use static_cast instead of dynamic_cast)
	//  since exact type can be determined
	SampleSet* tmp = NULL;
	data::Type outputType = mapping::extension[ext];

	if (outputType != rep->type()) {
		// output type differs from current representation type: conversion necessary
		switch (outputType) {
			case data::raw:
				switch (rep->type()) {
					case data::segmented:
						tmp = new RawSampleSet<cnvalue>(*static_cast<SegmentedSampleSet<cnvalue>*>(rep));
						swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case data::raw_ascn:
				switch (rep->type()) {
					case data::segmented_ascn:
						tmp = new RawSampleSet<alleles_cn>(*static_cast<SegmentedSampleSet<alleles_cn>*>(rep));
						swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case data::segmented:
				switch (rep->type()) {
					case data::raw:
						tmp = new SegmentedSampleSet<cnvalue>(*static_cast<RawSampleSet<cnvalue>*>(rep));
						swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case data::segmented_ascn:
				switch (rep->type()) {
					case data::raw_ascn:
						tmp = new SegmentedSampleSet<alleles_cn>(*static_cast<RawSampleSet<alleles_cn>*>(rep));
						swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			default:
				throw invalid_conversion("Cannot determine file type from file name extension.");
		}
	}
	
	rep->_write(file);
}
