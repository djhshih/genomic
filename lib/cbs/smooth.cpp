#include "cbs/smooth.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/math/distributions/normal.hpp>

namespace cbs {
namespace {

double inflfact(double trim) {
    using boost::math::normal_distribution;
    using boost::math::quantile;
    if (!(trim >= 0.0 && trim < 0.5)) {
        throw std::invalid_argument("trim must satisfy 0 <= trim < 0.5");
    }
    normal_distribution<double> nd;
    const double a = quantile(nd, 1.0 - trim);
    const int ngrid = 10000;
    const double step = (2.0 * a) / ngrid;
    double sum = 0.0;
    for (int i = 0; i < ngrid; ++i) {
        const double left = -a + i * step;
        const double right = left + step;
        const double x = 0.5 * (left + right);
        sum += x * x * pdf(nd, x) / (1.0 - 2.0 * trim);
    }
    return 1.0 / (sum * step);
}

double trimmed_variance(const std::vector<double>& genomdat, double trim) {
    const std::size_t n = genomdat.size();
    if (n < 2) return 0.0;
    const long long n_keep = std::llround((1.0 - 2.0 * trim) * static_cast<double>(n - 1));
    if (n_keep <= 0) return 0.0;
    std::vector<double> diffs;
    diffs.reserve(n - 1);
    for (std::size_t i = 1; i < n; ++i) diffs.push_back(std::abs(genomdat[i] - genomdat[i - 1]));
    std::sort(diffs.begin(), diffs.end());
    double ss = 0.0;
    for (long long i = 0; i < n_keep; ++i) ss += diffs[static_cast<std::size_t>(i)] * diffs[static_cast<std::size_t>(i)];
    return inflfact(trim) * (ss / (2.0 * static_cast<double>(n_keep)));
}

std::vector<int> finite_chrom_frequencies(const std::vector<int>& finite_chrom) {
    std::vector<int> cfrq;
    if (finite_chrom.empty()) return cfrq;
    int run = 1;
    for (std::size_t i = 1; i < finite_chrom.size(); ++i) {
        if (finite_chrom[i] == finite_chrom[i - 1]) {
            ++run;
        } else {
            cfrq.push_back(run);
            run = 1;
        }
    }
    cfrq.push_back(run);
    return cfrq;
}

double neighborhood_median(const std::vector<double>& gdat, int ilo, int ihi) {
    std::vector<double> xnbhd;
    xnbhd.reserve(static_cast<std::size_t>(ihi - ilo + 1));
    for (int j = ilo; j <= ihi; ++j) xnbhd.push_back(gdat[static_cast<std::size_t>(j)]);
    std::sort(xnbhd.begin(), xnbhd.end());
    const int k1 = static_cast<int>(xnbhd.size());
    const int j1 = k1 / 2;
    if (k1 == 2 * j1) {
        return (xnbhd[static_cast<std::size_t>(j1 - 1)] + xnbhd[static_cast<std::size_t>(j1)]) / 2.0;
    }
    return xnbhd[static_cast<std::size_t>(j1)];
}

std::vector<double> smooth_lr_kernel(const std::vector<double>& gdat,
                                     const std::vector<int>& cfrq,
                                     int k,
                                     double oSD,
                                     double sSD) {
    std::vector<double> sgdat(gdat.size());
    int cilo = 0;
    int cihi = -1;
    for (int freq : cfrq) {
        cihi += freq;
        for (int i = cilo; i <= cihi; ++i) {
            const int ilo = std::max(cilo, i - k);
            const int ihi = std::min(cihi, i + k);
            double mxnbd = 100.0 * oSD;
            double mnnbd = 100.0 * oSD;
            bool keep = false;
            for (int j = ilo; j <= ihi; ++j) {
                if (j == i) continue;
                const double distij = gdat[static_cast<std::size_t>(i)] - gdat[static_cast<std::size_t>(j)];
                if (std::abs(distij) <= oSD) {
                    sgdat[static_cast<std::size_t>(i)] = gdat[static_cast<std::size_t>(i)];
                    keep = true;
                    break;
                }
                if (distij < mxnbd) mxnbd = distij;
                if (-distij < mnnbd) mnnbd = -distij;
            }
            if (keep) continue;
            if ((mxnbd <= 0.0) && (mnnbd <= 0.0)) {
                sgdat[static_cast<std::size_t>(i)] = gdat[static_cast<std::size_t>(i)];
                continue;
            }
            const double xmed = neighborhood_median(gdat, ilo, ihi);
            if (mxnbd > 0.0) sgdat[static_cast<std::size_t>(i)] = xmed + sSD;
            if (mnnbd > 0.0) sgdat[static_cast<std::size_t>(i)] = xmed - sSD;
        }
        cilo += freq;
    }
    return sgdat;
}

} // namespace

std::vector<double> smooth(const std::vector<double>& values,
                               const std::vector<int>& chrom,
                               int smooth_region,
                               double outlier_sd_scale,
                               double smooth_sd_scale,
                               double trim) {
    if (values.size() != chrom.size()) throw std::invalid_argument("values and chrom must have same length");
    if (smooth_region < 0) throw std::invalid_argument("smooth_region must be non-negative");

    std::vector<double> out = values;
    std::vector<std::size_t> finite_idx;
    std::vector<double> finite_vals;
    std::vector<int> finite_chrom;
    finite_idx.reserve(values.size());
    finite_vals.reserve(values.size());
    finite_chrom.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (std::isfinite(values[i])) {
            finite_idx.push_back(i);
            finite_vals.push_back(values[i]);
            finite_chrom.push_back(chrom[i]);
        }
    }
    if (finite_vals.size() < 2) return out;

    const double tvar = trimmed_variance(finite_vals, trim);
    if (!(std::isfinite(tvar)) || tvar < 0.0) return out;
    const double trimmed_sd = std::sqrt(tvar);
    const double outlier_sd = outlier_sd_scale * trimmed_sd;
    const double smooth_sd = smooth_sd_scale * trimmed_sd;
    const std::vector<int> cfrq = finite_chrom_frequencies(finite_chrom);
    const std::vector<double> smoothed = smooth_lr_kernel(finite_vals, cfrq, smooth_region, outlier_sd, smooth_sd);
    for (std::size_t i = 0; i < finite_idx.size(); ++i) out[finite_idx[i]] = smoothed[i];
    return out;
}

std::vector<std::vector<double>> smooth_matrix(const std::vector<std::vector<double>>& samples,
                                                   const std::vector<int>& chrom,
                                                   int smooth_region,
                                                   double outlier_sd_scale,
                                                   double smooth_sd_scale,
                                                   double trim) {
    std::vector<std::vector<double>> out;
    out.reserve(samples.size());
    for (const auto& sample : samples) {
        out.push_back(smooth(sample, chrom, smooth_region, outlier_sd_scale, smooth_sd_scale, trim));
    }
    return out;
}

} // namespace cbs
