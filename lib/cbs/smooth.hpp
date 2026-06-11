#ifndef CNA_LIB_CBS_SMOOTH_HPP
#define CNA_LIB_CBS_SMOOTH_HPP

#include <vector>

namespace cbs {

std::vector<double> smooth(const std::vector<double>& values,
                               const std::vector<int>& chrom,
                               int smooth_region = 10,
                               double outlier_sd_scale = 4.0,
                               double smooth_sd_scale = 2.0,
                               double trim = 0.025);

std::vector<std::vector<double>> smooth_matrix(const std::vector<std::vector<double>>& samples,
                                                   const std::vector<int>& chrom,
                                                   int smooth_region = 10,
                                                   double outlier_sd_scale = 4.0,
                                                   double smooth_sd_scale = 2.0,
                                                   double trim = 0.025);

} // namespace cbs

#endif
