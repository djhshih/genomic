#include "GenericSampleSet.hpp"

void cna::GenericSampleSet::_read(std::fstream& file)
{
	const std::string& fileName = Base::fileName;
	cna::marker::Set* markers = Base::markers;
	
	std::string ext = cna::name::fileext(fileName);
	switch (cna::mapping::extension[ext]) {
		case cna::data::raw:
			rep = new cna::RawSampleSet<cnvalue>(markers);
			break;
		case cna::data::raw_ascn:
			rep = new cna::RawSampleSet<alleles_cn>(markers);
			break;
		case cna::data::segmented:
			rep = new cna::SegmentedSampleSet<cnvalue>(markers);
			break;
		case cna::data::segmented_ascn:
			rep = new cna::SegmentedSampleSet<alleles_cn>(markers);
			break;
		default:
			throw invalid_conversion("Cannot determine file type from file name extension.");
	}
	rep->_read(file);
}

void cna::GenericSampleSet::_write(std::fstream& file)
{
	const std::string& fileName = Base::fileName;
	
	std::string ext = cna::name::fileext(fileName);
	
	// cast $rep to appropriate type
	// runtime checking is skipped (i.e. use static_cast instead of dynamic_cast)
	//  since exact type can be determined
	cna::SampleSet* tmp = NULL;
	cna::data::Type outputType = cna::mapping::extension[ext];

	if (outputType != rep->type()) {
		// output type differs from current representation type: conversion necessary
		switch (outputType) {
			case cna::data::raw:
				switch (rep->type()) {
					case cna::data::segmented:
						tmp = new cna::RawSampleSet<cnvalue>(*static_cast<cna::SegmentedSampleSet<cnvalue>*>(rep));
						std::swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case cna::data::raw_ascn:
				switch (rep->type()) {
					case cna::data::segmented_ascn:
						tmp = new cna::RawSampleSet<alleles_cn>(*static_cast<cna::SegmentedSampleSet<alleles_cn>*>(rep));
						std::swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case cna::data::segmented:
				switch (rep->type()) {
					case cna::data::raw:
						tmp = new cna::SegmentedSampleSet<cnvalue>(*static_cast<cna::RawSampleSet<cnvalue>*>(rep));
						std::swap(rep, tmp);
						delete tmp;
						break;
					default:
						throw invalid_conversion("Conversion not supported.");
				}
				break;
			case cna::data::segmented_ascn:
				switch (rep->type()) {
					case cna::data::raw_ascn:
						tmp = new cna::SegmentedSampleSet<alleles_cn>(*static_cast<cna::RawSampleSet<alleles_cn>*>(rep));
						std::swap(rep, tmp);
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
