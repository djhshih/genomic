#ifndef CNA_LIB_CBS_CBS_HPP
#define CNA_LIB_CBS_CBS_HPP

#include <array>
#include <cstddef>
#include <random>
#include <vector>

namespace cbs {

struct ChangePointResult {
    int ncpt = 0;
    std::array<int, 2> icpt{{0, 0}};
    std::array<int, 2> iseg{{0, 0}};
    double ostat = 0.0;
};

struct BinarySegmentationResult {
    double statistic = 0.0;
    int start = 0;
    int end = 0;
};

struct SegmentationResult {
    std::vector<int> lengths;
    std::vector<double> means;
};

double tailp(double b, double delta, int m, int ngrid, double tol);
double btailp(double b, int m, int ng, double tol);

double btmax(const std::vector<double>& x);
BinarySegmentationResult tmaxo(const std::vector<double>& x, double tss, int al0, bool ibin);
double tmaxp(const std::vector<double>& px, double tss, int al0, bool ibin);
double htmaxp(const std::vector<double>& px, double tss, int k, int al0, bool ibin);

double tpermp(int n1, int n2, int n, const double* x, std::vector<double>& px, int nperm,
              std::mt19937_64& rng);
void xperm(const std::vector<double>& x, std::vector<double>& px, std::mt19937_64& rng);

void wxperm(const std::vector<double>& x,
            std::vector<double>& px,
            const std::vector<double>& rwts,
            std::mt19937_64& rng);
double wtpermp(int n1,
               int n2,
               int n,
               const double* x,
               std::vector<double>& px,
               const double* wts,
               const double* rwts,
               int nperm,
               std::mt19937_64& rng);
void getmncwt(int n, const std::vector<double>& cwts, int k, std::vector<double>& mncwt, double& delta);
BinarySegmentationResult wtmaxo(const std::vector<double>& x,
                                const std::vector<double>& wts,
                                double tss,
                                const std::vector<double>& cwts,
                                int al0);
double wtmaxp(const std::vector<double>& px,
              const std::vector<double>& wts,
              const std::vector<double>& cwts,
              int al0);
double hwtmaxp(const std::vector<double>& px,
               const std::vector<double>& wts,
               const std::vector<double>& cwts,
               const std::vector<double>& mncwt,
               int k,
               int al0);

ChangePointResult fndcpt(const std::vector<double>& x,
                         double tss,
                         int nperm,
                         double cpval,
                         bool ibin,
                         bool hybrid,
                         int al0,
                         int hk,
                         double delta,
                         int ngrid,
                         const std::vector<int>& sbdry,
                         double tol,
                         std::mt19937_64& rng);
ChangePointResult wfindcpt(const std::vector<double>& x,
                           double tss,
                           const std::vector<double>& wts,
                           const std::vector<double>& rwts,
                           const std::vector<double>& cwts,
                           int nperm,
                           double cpval,
                           bool hybrid,
                           int al0,
                           int hk,
                           double delta,
                           int ngrid,
                           const std::vector<int>& sbdry,
                           double tol,
                           std::mt19937_64& rng);

SegmentationResult segment(const std::vector<double>& x,
                           bool ibin,
                           double alpha,
                           int nperm,
                           bool hybrid,
                           int min_width,
                           int kmax,
                           int nmin,
                           double eta,
                           const std::vector<int>& sbdry,
                           double tol,
                           std::mt19937_64& rng,
                           bool undo_prune = false,
                           double undo_prune_cutoff = 0.05);

SegmentationResult segment_weighted(const std::vector<double>& x,
                                    const std::vector<double>& weights,
                                    double alpha,
                                    int nperm,
                                    bool hybrid,
                                    int min_width,
                                    int kmax,
                                    int nmin,
                                    double eta,
                                    const std::vector<int>& sbdry,
                                    double tol,
                                    std::mt19937_64& rng,
                                    bool undo_prune = false,
                                    double undo_prune_cutoff = 0.05);

} // namespace cbs

#endif
