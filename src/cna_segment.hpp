#ifndef cna_segment_h
#define cna_segment_h

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
#include "global.hpp"
#include "cna_common.hpp"
#include "SampleSets.hpp"
#include "cbs/smooth.hpp"
#include "cbs/CBS.hpp"

class Segment : public Command {
public:
	Segment()
	: Command("smooth and CBS-segment a raw copy-number matrix") {
		opts.add_options()
			("help", "print help message")
			("input,i", po::value<std::string>(), "raw sample matrix file")
			("output,o", po::value<std::string>(), "output segmentation file")
			("format,f", po::value<std::string>(), "input file format [default: determined from file extension]")
			("alpha", po::value<double>(), "CBS alpha [default: 0.01]")
			("nperm", po::value<int>(), "CBS permutations [default: 200]")
			("min_width", po::value<int>(), "CBS minimum segment width [default: 2]")
			("kmax", po::value<int>(), "CBS kmax [default: 25]")
			("nmin", po::value<int>(), "CBS nmin [default: 200]")
			("eta", po::value<double>(), "CBS eta [default: 0.05]")
			("trim", po::value<double>(), "smooth trim [default: 0.025]")
			("smooth_region", po::value<int>(), "smooth neighborhood radius [default: 10]")
			("outlier_sd_scale", po::value<double>(), "smooth outlier SD scale [default: 4.0]")
			("smooth_sd_scale", po::value<double>(), "smooth replacement SD scale [default: 2.0]")
			("hybrid", po::value<bool>(), "use hybrid CBS p-values [default: false]")
			("undo_prune", po::value<bool>(), "apply prune undo [default: false]")
			("undo_prune_cutoff", po::value<double>(), "prune cutoff [default: 0.05]")
			;
		popts.add("input", 1).add("output", 1);
	}

	void run() override {
		if (vm.count("help")) {
			std::cout << "usage:  " << progname << " segment [options] <raw sample matrix file> <output segmentation file>" << std::endl;
			std::cout << opts << std::endl;
			return;
		}

		getOptions();
		if (inputType != cna::data::raw) {
			throw std::invalid_argument("segment command currently supports raw log-ratio matrices only.");
		}

		cna::RawSampleSet<rvalue> raw;
		raw.read(inputFileName);
		ensure_log_scale(raw);
		cna::SegmentedSampleSet<rvalue> segmented = segment_raw(raw);
		segmented.write(outputFileName);
	}

private:
	std::string inputFileName, outputFileName;
	cna::data::Type inputType = cna::data::invalid;
	double alpha = 0.01;
	int nperm = 200;
	int minWidth = 2;
	int kmax = 25;
	int nmin = 200;
	double eta = 0.05;
	double trim = 0.025;
	int smoothRegion = 10;
	double outlierSdScale = 4.0;
	double smoothSdScale = 2.0;
	bool hybrid = false;
	bool undoPrune = false;
	double undoPruneCutoff = 0.05;

	void getOptions() {
		if (vm.count("input")) inputFileName = vm["input"].as<std::string>();
		else throw std::invalid_argument("Input file not specified.");

		if (vm.count("format")) inputType = cna::mapping::extension[vm["format"].as<std::string>()];
		else inputType = cna::mapping::extension[cna::name::fileext(inputFileName)];
		if (inputType == cna::data::invalid) {
			throw std::invalid_argument("Invalid input format type for input file '" + inputFileName + "'.");
		}

		if (vm.count("output")) outputFileName = vm["output"].as<std::string>();
		else outputFileName = cna::name::filestem(inputFileName) + ".seg";

		if (vm.count("alpha")) alpha = vm["alpha"].as<double>();
		if (vm.count("nperm")) nperm = vm["nperm"].as<int>();
		if (vm.count("min_width")) minWidth = vm["min_width"].as<int>();
		if (vm.count("kmax")) kmax = vm["kmax"].as<int>();
		if (vm.count("nmin")) nmin = vm["nmin"].as<int>();
		if (vm.count("eta")) eta = vm["eta"].as<double>();
		if (vm.count("trim")) trim = vm["trim"].as<double>();
		if (vm.count("smooth_region")) smoothRegion = vm["smooth_region"].as<int>();
		if (vm.count("outlier_sd_scale")) outlierSdScale = vm["outlier_sd_scale"].as<double>();
		if (vm.count("smooth_sd_scale")) smoothSdScale = vm["smooth_sd_scale"].as<double>();
		if (vm.count("hybrid")) hybrid = vm["hybrid"].as<bool>();
		if (vm.count("undo_prune")) undoPrune = vm["undo_prune"].as<bool>();
		if (vm.count("undo_prune_cutoff")) undoPruneCutoff = vm["undo_prune_cutoff"].as<double>();
	}

	static void ensure_log_scale(cna::RawSampleSet<rvalue>& raw) {
		bool has_neg = false, has_pos = false;
		for (const auto* sample_ptr : raw.getSamples()) {
			auto* sample = const_cast<cna::RawSampleSet<rvalue>::RawSample*>(sample_ptr);
			for (auto chr_it = sample->begin(); chr_it != sample->end(); ++chr_it) {
				for (auto v_it = chr_it->begin(); v_it != chr_it->end(); ++v_it) {
					const double v = *v_it;
					if (!std::isfinite(v)) continue;
					if (v < 0) has_neg = true;
					if (v > 0) has_pos = true;
				}
			}
		}
		if (!(has_neg && has_pos)) {
			throw std::invalid_argument("Input does not appear to be in log scale: expected both negative and positive values.");
		}
	}

	cna::SegmentedSampleSet<rvalue> segment_raw(cna::RawSampleSet<rvalue>& raw) const {
		cna::SegmentedSampleSet<rvalue> out(raw.marker_set());
		std::mt19937_64 rng(1);
		std::vector<int> sbdry((nperm + 1) * (nperm + 2) / 2 + 2, nperm + 1);

		for (const auto* sample_ptr : raw.getSamples()) {
			auto* sample_it = const_cast<cna::RawSampleSet<rvalue>::RawSample*>(sample_ptr);
			auto* out_sample = out.create(sample_it->name);
			for (std::size_t chri = 0; chri < sample_it->size(); ++chri) {
				auto& chr = (*sample_it)[static_cast<chromid>(chri)];
				std::vector<double> x(chr.begin(), chr.end());
				if (x.empty()) continue;
				std::vector<int> chrom(x.size(), static_cast<int>(chri + 1));
				const std::vector<double> smoothed = cbs::smooth(x, chrom, smoothRegion, outlierSdScale, smoothSdScale, trim);
				const auto seg = cbs::segment(smoothed, false, alpha, nperm, hybrid, minWidth, kmax, nmin, eta, sbdry, 1e-6, rng, undoPrune, undoPruneCutoff);

				std::size_t start_index = 0;
				for (std::size_t i = 0; i < seg.lengths.size(); ++i) {
					const std::size_t len = static_cast<std::size_t>(seg.lengths[i]);
					if (len == 0) continue;
					const std::size_t end_index = start_index + len - 1;
					cna::Segment<rvalue> s(static_cast<chromid>(chri + 1),
						raw.marker_set()->at(chri)[start_index]->pos,
						raw.marker_set()->at(chri)[end_index]->pos,
						static_cast<unsigned long>(len),
						seg.means[i]);
					out_sample->chromosome(static_cast<chromid>(chri))->push_back(s);
					start_index += len;
				}
			}
		}
		return out;
	}
};

#endif
