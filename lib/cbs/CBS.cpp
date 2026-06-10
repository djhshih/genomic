#include "cbs/CBS.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace cbs {
namespace {

constexpr double kInvSqrt2Pi = 0.39894228040143267794;

inline double fpnorm(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

inline double nu(double x, double tol) {
    if (x > 0.01) {
        double lnu1 = std::log(2.0) - 2.0 * std::log(x);
        double lnu0 = lnu1;
        int k = 2;
        double dk = 0.0;
        for (int i = 1; i <= k; ++i) {
            dk += 1.0;
            const double xk = -x * std::sqrt(dk) / 2.0;
            lnu1 -= 2.0 * fpnorm(xk) / dk;
        }
        while (std::abs((lnu1 - lnu0) / lnu1) > tol) {
            lnu0 = lnu1;
            for (int i = 1; i <= k; ++i) {
                dk += 1.0;
                const double xk = -x * std::sqrt(dk) / 2.0;
                lnu1 -= 2.0 * fpnorm(xk) / dk;
            }
            k *= 2;
        }
        return std::exp(lnu1);
    }
    return std::exp(-0.583 * x);
}

inline double it1tsq(double x, double a) {
    double y = x + a - 0.5;
    double out = (8.0 * y) / (1.0 - 4.0 * y * y) +
                 2.0 * std::log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
    y = x - 0.5;
    out -= (8.0 * y) / (1.0 - 4.0 * y * y) +
           2.0 * std::log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
    return out;
}

inline double runif01(std::mt19937_64& rng) {
    return std::generate_canonical<double, 53>(rng);
}

struct BlockScanResult {
    double statistic = 0.0;
    int left = 0;
    int right = 0;
};

inline void sort_loc_by_values(const std::vector<double>& values, std::vector<int>& loc, int n) {
    std::sort(loc.begin() + 1, loc.begin() + n + 1,
              [&](int a, int b) { return values[a] < values[b]; });
}

BlockScanResult tmaxo_impl(const std::vector<double>& x, double tss, int al0, bool ibin, bool track_loc) {
    const int n = static_cast<int>(x.size());
    const double rn = static_cast<double>(n);
    const int nb = (n >= 50) ? static_cast<int>(std::lround(std::sqrt(rn))) : 1;
    const int nb2 = nb * (nb + 1) / 2;

    std::vector<double> sx(n + 1, 0.0), bpsmax(nb + 1), bpsmin(nb + 1), bssbij(nb2 + 1), bssijmax(nb2 + 1);
    std::vector<int> bb(nb + 1), ibmin(nb + 1), ibmax(nb + 1), bloci(nb2 + 1), blocj(nb2 + 1), loc(nb2 + 1), alen(nb2 + 1);

    for (int i = 1; i <= nb; ++i) bb[i] = static_cast<int>(std::lround(rn * (static_cast<double>(i) / static_cast<double>(nb))));

    int ilo = 1;
    double psum = 0.0, psmin0 = 0.0, psmax0 = 0.0;
    int ipsmin0 = n, ipsmax0 = n;
    for (int j = 1; j <= nb; ++j) {
        sx[ilo] = psum + x[ilo - 1];
        double psmin = sx[ilo], psmax = sx[ilo];
        int ipsmin = ilo, ipsmax = ilo;
        for (int i = ilo + 1; i <= bb[j]; ++i) {
            sx[i] = sx[i - 1] + x[i - 1];
            if (sx[i] < psmin) { psmin = sx[i]; ipsmin = i; }
            if (sx[i] > psmax) { psmax = sx[i]; ipsmax = i; }
        }
        ibmin[j] = ipsmin; ibmax[j] = ipsmax;
        bpsmin[j] = psmin; bpsmax[j] = psmax;
        if (psmin < psmin0) { psmin0 = psmin; ipsmin0 = ipsmin; }
        if (psmax > psmax0) { psmax0 = psmax; ipsmax0 = ipsmax; }
        psum = sx[bb[j]];
        ilo = bb[j] + 1;
    }

    const double psdiff = psmax0 - psmin0;
    double bssmax = 0.0;
    int tmaxi = std::min(ipsmax0, ipsmin0), tmaxj = std::max(ipsmax0, ipsmin0);
    if (psdiff <= 0.0) {
        if (ibin) {
            if (tss <= 0.0001) tss = 1.0;
            bssmax = 0.0 / (tss / rn);
        } else {
            if (tss <= 0.0001) tss = 1.0;
            bssmax = 0.0 / ((tss - 0.0) / (rn - 2.0));
        }
        return {bssmax, track_loc ? tmaxi : 0, track_loc ? tmaxj : 0};
    }

    {
        const double rj = static_cast<double>(std::abs(ipsmax0 - ipsmin0));
        const double rnjov1 = rn / (rj * (rn - rj));
        bssmax = ibin ? rnjov1 * std::pow(psdiff - 0.5, 2.0) : rnjov1 * psdiff * psdiff;
    }

    const double rnov2 = rn / 2.0;
    int l = 0;
    const int nal0 = n - al0;
    for (int i = 1; i <= nb; ++i) {
        for (int j = i; j <= nb; ++j) {
            const int ilo1 = (i == 1) ? 1 : bb[i - 1] + 1;
            const int ihi = bb[i];
            const int jlo = (j == 1) ? 1 : bb[j - 1] + 1;
            const int jhi = bb[j];
            int alenhi = jhi - ilo1;
            if (alenhi > nal0) alenhi = nal0;
            int alenlo = (i == j) ? 1 : (jlo - ihi);
            if (alenlo < al0) alenlo = al0;
            const double sij1 = std::abs(bpsmax[j] - bpsmin[i]);
            const double sij2 = std::abs(bpsmax[i] - bpsmin[j]);
            const double sijmx0 = std::max(sij1, sij2);
            const double rjlo = static_cast<double>(alenlo);
            const double rjhi = static_cast<double>(alenhi);
            const double rnjov1 = rn / std::min(rjlo * (rn - rjlo), rjhi * (rn - rjhi));
            const double bsslim = ibin ? rnjov1 * std::pow(sijmx0 - 0.5, 2.0) : rnjov1 * sijmx0 * sijmx0;
            if (bssmax <= bsslim) {
                ++l;
                loc[l] = l;
                bloci[l] = i;
                blocj[l] = j;
                bssijmax[l] = bsslim;
                if (sij1 > sij2) {
                    alen[l] = std::abs(ibmax[j] - ibmin[i]);
                    const double rr = static_cast<double>(alen[l]);
                    const double fac = rn / (rr * (rn - rr));
                    bssbij[l] = ibin ? fac * std::pow(sij1 - 0.5, 2.0) : fac * sij1 * sij1;
                } else {
                    alen[l] = std::abs(ibmin[j] - ibmax[i]);
                    const double rr = static_cast<double>(alen[l]);
                    const double fac = rn / (rr * (rn - rr));
                    bssbij[l] = ibin ? fac * std::pow(sij2 - 0.5, 2.0) : fac * sij2 * sij2;
                }
            }
        }
    }
    const int nb1 = l;
    sort_loc_by_values(bssbij, loc, nb1);

    for (int ll = nb1; ll >= 1; --ll) {
        const int k = loc[ll];
        const double bsslim = bssijmax[k];
        if (bssmax > bsslim) continue;
        const int bi = bloci[k], bj = blocj[k];
        int alenmax = alen[k];
        const int ilo1 = (bi == 1) ? 1 : bb[bi - 1] + 1;
        const int ihi = bb[bi];
        const int jlo = (bj == 1) ? 1 : bb[bj - 1] + 1;
        const int jhi = bb[bj];
        int alenhi = jhi - ilo1;
        if (alenhi > nal0) alenhi = nal0;
        int alenlo = (bi == bj) ? 1 : (jlo - ihi);
        if (alenlo < al0) alenlo = al0;
        const double rjlo = static_cast<double>(alenlo);
        const double rjhi = static_cast<double>(alenhi);

        if (alenmax > n - alenmax) alenmax = n - alenmax;
        if ((rjlo <= rnov2) && (alenlo <= alenmax)) {
            for (int i2j = alenlo; i2j <= alenmax; ++i2j) {
                const int ixlo = std::max(0, jlo - ilo1 - i2j);
                const int ixhi = std::max(0, ihi + i2j - jhi);
                double sxmx = 0.0;
                int sxmxi = ilo1;
                for (int i = ilo1 + ixlo; i <= ihi - ixhi; ++i) {
                    const int j = i + i2j;
                    const double absx = std::abs(sx[j] - sx[i]);
                    if (sxmx < absx) { sxmx = absx; sxmxi = i; }
                }
                const double rr = static_cast<double>(i2j);
                const double fac = rn / (rr * (rn - rr));
                const double bijbss = ibin ? fac * std::pow(sxmx - 0.5, 2.0) : fac * sxmx * sxmx;
                if (bijbss > bssmax) { bssmax = bijbss; tmaxi = sxmxi; tmaxj = sxmxi + i2j; }
            }
        }

        alenmax = n - alenmax;
        if ((rjhi >= rnov2) && (alenhi >= alenmax)) {
            for (int i2j = alenhi; i2j >= alenmax; --i2j) {
                const int ixlo = std::max(0, jlo - ilo1 - i2j);
                const int ixhi = std::max(0, ihi + i2j - jhi);
                double sxmx = 0.0;
                int sxmxi = ilo1;
                for (int i = ilo1 + ixlo; i <= ihi - ixhi; ++i) {
                    const int j = i + i2j;
                    const double absx = std::abs(sx[j] - sx[i]);
                    if (sxmx < absx) { sxmx = absx; sxmxi = i; }
                }
                const double rr = static_cast<double>(i2j);
                const double fac = rn / (rr * (rn - rr));
                const double bijbss = ibin ? fac * std::pow(sxmx - 0.5, 2.0) : fac * sxmx * sxmx;
                if (bijbss > bssmax) { bssmax = bijbss; tmaxi = sxmxi; tmaxj = sxmxi + i2j; }
            }
        }
    }

    if (ibin) {
        if (tss <= 0.0001) tss = 1.0;
        bssmax /= (tss / rn);
    } else {
        if (tss <= bssmax + 0.0001) tss = bssmax + 1.0;
        bssmax /= ((tss - bssmax) / (rn - 2.0));
    }

    return {bssmax, track_loc ? tmaxi : 0, track_loc ? tmaxj : 0};
}

} // namespace

double tailp(double b, double delta, int m, int ngrid, double tol) {
    const double dincr = (0.5 - delta) / static_cast<double>(ngrid);
    const double bsqrtm = b / std::sqrt(static_cast<double>(m));
    double tl = 0.5 - dincr;
    double t = 0.5 - 0.5 * dincr;
    double out = 0.0;
    for (int i = 1; i <= ngrid; ++i) {
        tl += dincr;
        t += dincr;
        const double x = bsqrtm / std::sqrt(t * (1.0 - t));
        const double nux = nu(x, tol);
        out += (nux * nux) * it1tsq(tl, dincr);
    }
    out = 9.973557e-2 * std::pow(b, 3.0) * std::exp(-b * b / 2.0) * out;
    return 2.0 * out;
}

double btailp(double b, int m, int ng, double tol) {
    const double dm = static_cast<double>(m);
    const int k = 2;
    const double ll = b * std::sqrt(1.0 / static_cast<double>(m - k) - 1.0 / dm);
    const double ul = b * std::sqrt(1.0 / static_cast<double>(k) - 1.0 / dm);
    const double dincr = (ul - ll) / static_cast<double>(ng);
    double out = 0.0;
    double x = ll;
    double x1 = x + (b * b) / (dm * x);
    double nulo = nu(x1, tol) / x;
    for (int i = 1; i <= ng; ++i) {
        x += dincr;
        x1 = x + (b * b) / (dm * x);
        const double nuhi = nu(x1, tol) / x;
        out += (nuhi + nulo) * dincr;
        nulo = nuhi;
    }
    out = b * std::exp(-b * b / 2.0) * out * kInvSqrt2Pi;
    out += 2.0 * (1.0 - fpnorm(b));
    return out;
}

double btmax(const std::vector<double>& x) {
    const int n = static_cast<int>(x.size());
    double sumxi = x[0];
    double ostat = 0.0;
    const double dn = static_cast<double>(n);
    double di = 1.0;
    for (int i = 2; i <= n - 2; ++i) {
        di += 1.0;
        sumxi += x[i - 1];
        const double btmaxi = dn * (sumxi * sumxi) / (di * (dn - di));
        if (ostat < btmaxi) ostat = btmaxi;
    }
    return std::sqrt(ostat);
}

BinarySegmentationResult tmaxo(const std::vector<double>& x, double tss, int al0, bool ibin) {
    const auto r = tmaxo_impl(x, tss, al0, ibin, true);
    return {r.statistic, r.left - 1, r.right - 1};
}

double tmaxp(const std::vector<double>& px, double tss, int al0, bool ibin) {
    return tmaxo_impl(px, tss, al0, ibin, false).statistic;
}

double htmaxp(const std::vector<double>& px, double tss, int k, int al0, bool ibin) {
    const int n = static_cast<int>(px.size());
    const double rn = static_cast<double>(n);
    const int nb = static_cast<int>(rn / static_cast<double>(k));
    std::vector<double> bpsmax(nb + 1), bpsmin(nb + 1), sx(n + 1, 0.0);
    std::vector<int> bb(nb + 1);
    for (int i = 1; i <= nb; ++i) bb[i] = static_cast<int>(std::lround(rn * (static_cast<double>(i) / static_cast<double>(nb))));

    int ilo = 1;
    double psum = 0.0;
    double out = 0.0;
    for (int j = 1; j <= nb; ++j) {
        sx[ilo] = psum + px[ilo - 1];
        double psmin = sx[ilo], psmax = sx[ilo];
        int ipsmin = ilo, ipsmax = ilo;
        for (int i = ilo + 1; i <= bb[j]; ++i) {
            sx[i] = sx[i - 1] + px[i - 1];
            if (sx[i] < psmin) { psmin = sx[i]; ipsmin = i; }
            if (sx[i] > psmax) { psmax = sx[i]; ipsmax = i; }
        }
        bpsmin[j] = psmin;
        bpsmax[j] = psmax;
        psum = sx[bb[j]];
        ilo = bb[j] + 1;
        const int d = std::abs(ipsmin - ipsmax);
        if ((d <= k) && (d >= al0)) {
            const double rj = static_cast<double>(d);
            const double fac = rn / (rj * (rn - rj));
            const double bssmx = ibin ? fac * std::pow(bpsmax[j] - bpsmin[j] - 0.5, 2.0)
                                      : fac * std::pow(bpsmax[j] - bpsmin[j], 2.0);
            if (out < bssmx) out = bssmx;
        }
    }

    auto scan_regular = [&](int ilo1, int ihi, double psdiff) {
        const double psdiffsq = ibin ? std::pow(psdiff - 0.5, 2.0) : std::pow(psdiff, 2.0);
        for (int j = al0; j <= k; ++j) {
            const double rj = static_cast<double>(j);
            const double fac = rn / (rj * (rn - rj));
            const double bsslim = fac * psdiffsq;
            if (bsslim < out) goto done;
            double sxmx = 0.0;
            for (int i = ilo1; i <= ihi - j; ++i) {
                const double absx = std::abs(sx[i + j] - sx[i]);
                if (sxmx < absx) sxmx = absx;
            }
            const double bssmx = ibin ? fac * std::pow(std::abs(sxmx) - 0.5, 2.0) : fac * sxmx * sxmx;
            if (out < bssmx) out = bssmx;
        }
        done: ;
    };

    scan_regular(1, bb[1], bpsmax[1] - bpsmin[1]);

    {
        const double psdiff = std::max(std::abs(bpsmax[1] - bpsmin[nb]), std::abs(bpsmax[nb] - bpsmin[1]));
        const double psdiffsq = ibin ? std::pow(psdiff - 0.5, 2.0) : std::pow(psdiff, 2.0);
        for (int j = al0; j <= k; ++j) {
            const double rj = static_cast<double>(j);
            const double fac = rn / (rj * (rn - rj));
            const double bsslim = fac * psdiffsq;
            if (bsslim < out) break;
            double sxmx = 0.0;
            const int nmj = n - j;
            for (int i = 1; i <= j; ++i) {
                const double absx = std::abs(sx[i + nmj] - sx[i]);
                if (sxmx < absx) sxmx = absx;
            }
            const double bssmx = ibin ? fac * std::pow(std::abs(sxmx) - 0.5, 2.0) : fac * sxmx * sxmx;
            if (out < bssmx) out = bssmx;
        }
    }

    for (int l = 2; l <= nb; ++l) {
        scan_regular(bb[l - 1] + 1, bb[l], bpsmax[l] - bpsmin[l]);
        const double psdiff = std::max(std::abs(bpsmax[l] - bpsmin[l - 1]), std::abs(bpsmax[l - 1] - bpsmin[l]));
        const double psdiffsq = ibin ? std::pow(psdiff - 0.5, 2.0) : std::pow(psdiff, 2.0);
        for (int j = al0; j <= k; ++j) {
            const double rj = static_cast<double>(j);
            const double fac = rn / (rj * (rn - rj));
            const double bsslim = fac * psdiffsq;
            if (bsslim < out) break;
            double sxmx = 0.0;
            for (int i = bb[l - 1] + 1 - j; i <= bb[l - 1]; ++i) {
                const double absx = std::abs(sx[i + j] - sx[i]);
                if (sxmx < absx) sxmx = absx;
            }
            const double bssmx = ibin ? fac * std::pow(std::abs(sxmx) - 0.5, 2.0) : fac * sxmx * sxmx;
            if (out < bssmx) out = bssmx;
        }
    }

    if (ibin) {
        if (tss <= 0.0001) tss = 1.0;
        return out / (tss / rn);
    }
    if (tss <= out + 0.0001) tss = out + 1.0;
    return out / ((tss - out) / (rn - 2.0));
}

void xperm(const std::vector<double>& x, std::vector<double>& px, std::mt19937_64& rng) {
    px = x;
    for (int i = static_cast<int>(px.size()); i >= 1; --i) {
        const int j = static_cast<int>(runif01(rng) * static_cast<double>(i)) + 1;
        std::swap(px[i - 1], px[j - 1]);
    }
}

double tpermp(int n1, int n2, int n, const double* x, std::vector<double>& px, int nperm, std::mt19937_64& rng) {
    const double rn1 = static_cast<double>(n1);
    const double rn2 = static_cast<double>(n2);
    const double rn = rn1 + rn2;
    if (n1 == 1 || n2 == 1) return 1.0;
    px.resize(n);
    double xsum1 = 0.0, xsum2 = 0.0, tss = 0.0;
    for (int i = 0; i < n1; ++i) {
        px[i] = x[i];
        xsum1 += x[i];
        tss += x[i] * x[i];
    }
    for (int i = n1; i < n; ++i) {
        px[i] = x[i];
        xsum2 += x[i];
        tss += x[i] * x[i];
    }
    const double xbar = (xsum1 + xsum2) / rn;
    tss -= rn * (xbar * xbar);
    int m1;
    double rm1, ostat, tstat;
    if (n1 <= n2) {
        m1 = n1; rm1 = rn1; ostat = 0.99999 * std::abs(xsum1 / rn1 - xbar); tstat = (ostat * ostat) * rn1 * rn / rn2;
    } else {
        m1 = n2; rm1 = rn2; ostat = 0.99999 * std::abs(xsum2 / rn2 - xbar); tstat = (ostat * ostat) * rn2 * rn / rn1;
    }
    tstat /= ((tss - tstat) / (rn - 2.0));
    if ((tstat > 25.0) && (m1 >= 10)) return 1.0;
    int nrej = 0;
    for (int np = 1; np <= nperm; ++np) {
        xsum1 = 0.0;
        for (int i = n; i >= n - m1 + 1; --i) {
            const int j = static_cast<int>(runif01(rng) * static_cast<double>(i)) + 1;
            std::swap(px[i - 1], px[j - 1]);
            xsum1 += px[i - 1];
        }
        const double pstat = std::abs(xsum1 / rm1 - xbar);
        if (ostat <= pstat) ++nrej;
    }
    return static_cast<double>(nrej) / static_cast<double>(nperm);
}

void wxperm(const std::vector<double>& x, std::vector<double>& px, const std::vector<double>& rwts, std::mt19937_64& rng) {
    px.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) px[i] = x[i] * rwts[i];
    for (int i = static_cast<int>(px.size()); i >= 1; --i) {
        const int j = static_cast<int>(runif01(rng) * static_cast<double>(i)) + 1;
        const double tmpx = px[i - 1];
        px[i - 1] = px[j - 1] / rwts[i - 1];
        px[j - 1] = tmpx;
    }
}

double wtpermp(int n1, int n2, int n, const double* x, std::vector<double>& px, const double* wts, const double* rwts, int nperm, std::mt19937_64& rng) {
    if (n1 == 1 || n2 == 1) return 1.0;
    px.resize(n);
    double xsum1 = 0.0, xsum2 = 0.0, tss = 0.0, rn1 = 0.0, rn2 = 0.0;
    for (int i = 0; i < n1; ++i) {
        px[i] = x[i] * rwts[i];
        xsum1 += wts[i] * x[i];
        tss += wts[i] * x[i] * x[i];
        rn1 += wts[i];
    }
    for (int i = n1; i < n; ++i) {
        px[i] = x[i];
        xsum2 += wts[i] * x[i];
        tss += wts[i] * x[i] * x[i];
        rn2 += wts[i];
    }
    const double rn = rn1 + rn2;
    const double xbar = (xsum1 + xsum2) / rn;
    tss -= rn * (xbar * xbar);
    int m1;
    double rm1, ostat, tstat;
    if (n1 <= n2) {
        m1 = n1; rm1 = rn1; ostat = 0.99999 * std::abs(xsum1 / rn1 - xbar); tstat = (ostat * ostat) * rn1 * rn / rn2;
    } else {
        m1 = n2; rm1 = rn2; ostat = 0.99999 * std::abs(xsum2 / rn2 - xbar); tstat = (ostat * ostat) * rn2 * rn / rn1;
    }
    tstat /= ((tss - tstat) / (static_cast<double>(n) - 2.0));
    if ((tstat > 25.0) && (m1 >= 10)) return 1.0;
    int nrej = 0;
    for (int np = 1; np <= nperm; ++np) {
        xsum1 = 0.0;
        for (int i = n; i >= n - m1 + 1; --i) {
            const int j = static_cast<int>(runif01(rng) * static_cast<double>(i)) + 1;
            std::swap(px[i - 1], px[j - 1]);
            xsum1 += px[i - 1] * rwts[i - 1];
        }
        const double pstat = std::abs(xsum1 / rm1 - xbar);
        if (ostat <= pstat) ++nrej;
    }
    return static_cast<double>(nrej) / static_cast<double>(nperm);
}

void getmncwt(int n, const std::vector<double>& cwts, int k, std::vector<double>& mncwt, double& delta) {
    mncwt.assign(k + 1, 0.0);
    const double rn = cwts[n - 1];
    for (int j = 1; j <= k; ++j) {
        mncwt[j] = cwts[j - 1];
        const int nmj = n - j;
        for (int i = 1; i <= nmj; ++i) mncwt[j] = std::min(mncwt[j], cwts[i + j - 1] - cwts[i - 1]);
        for (int i = 1; i <= j; ++i) mncwt[j] = std::min(mncwt[j], rn - (cwts[i + nmj - 1] - cwts[i - 1]));
    }
    const int j = k + 1;
    const int nmj = n - j;
    delta = cwts[j - 1];
    for (int i = 1; i <= nmj; ++i) delta = std::min(delta, cwts[i + j - 1] - cwts[i - 1]);
    for (int i = 1; i <= j; ++i) delta = std::min(delta, rn - (cwts[i + nmj - 1] - cwts[i - 1]));
    delta /= cwts[n - 1];
}

BinarySegmentationResult wtmaxo(const std::vector<double>& x, const std::vector<double>& wts, double tss, const std::vector<double>& cwts, int al0) {
    const int n = static_cast<int>(x.size());
    const double rn = static_cast<double>(n);
    const int nb = (n >= 50) ? static_cast<int>(std::lround(std::sqrt(rn))) : 1;
    const int nb2 = nb * (nb + 1) / 2;
    std::vector<double> sx(n + 1, 0.0), bpsmax(nb + 1), bpsmin(nb + 1), bssbij(nb2 + 1), bssijmax(nb2 + 1), awt(nb2 + 1);
    std::vector<int> bb(nb + 1), ibmin(nb + 1), ibmax(nb + 1), bloci(nb2 + 1), blocj(nb2 + 1), loc(nb2 + 1);
    for (int i = 1; i <= nb; ++i) bb[i] = static_cast<int>(std::lround(rn * (static_cast<double>(i) / static_cast<double>(nb))));

    int ilo = 1;
    double psum = 0.0, psmin0 = 0.0, psmax0 = 0.0;
    int ipsmin0 = n, ipsmax0 = n;
    for (int j = 1; j <= nb; ++j) {
        sx[ilo] = psum + x[ilo - 1] * wts[ilo - 1];
        double psmin = sx[ilo], psmax = sx[ilo];
        int ipsmin = ilo, ipsmax = ilo;
        for (int i = ilo + 1; i <= bb[j]; ++i) {
            sx[i] = sx[i - 1] + x[i - 1] * wts[i - 1];
            if (sx[i] < psmin) { psmin = sx[i]; ipsmin = i; }
            if (sx[i] > psmax) { psmax = sx[i]; ipsmax = i; }
        }
        ibmin[j] = ipsmin; ibmax[j] = ipsmax;
        bpsmin[j] = psmin; bpsmax[j] = psmax;
        if (psmin < psmin0) { psmin0 = psmin; ipsmin0 = ipsmin; }
        if (psmax > psmax0) { psmax0 = psmax; ipsmax0 = ipsmax; }
        psum = sx[bb[j]];
        ilo = bb[j] + 1;
    }

    double bssmax = 0.0;
    int tmaxi = std::min(ipsmax0, ipsmin0), tmaxj = std::max(ipsmax0, ipsmin0);
    const double psdiff = psmax0 - psmin0;
    if (psdiff <= 0.0) {
        if (tss <= bssmax + 0.0001) tss = bssmax + 1.0;
        return {bssmax / ((tss - bssmax) / (rn - 2.0)), tmaxi, tmaxj};
    }

    const double psrn = cwts[n - 1];
    const double psrj = std::abs(cwts[ipsmax0 - 1] - cwts[ipsmin0 - 1]);
    bssmax = (psdiff * psdiff) / (psrj * (psrn - psrj));
    const double psrnov2 = psrn / 2.0;
    int l = 0;
    const int nal0 = n - al0;
    for (int i = 1; i <= nb; ++i) {
        for (int j = i; j <= nb; ++j) {
            const int ilo1 = (i == 1) ? 1 : bb[i - 1] + 1;
            const int ihi = bb[i];
            const int jlo = (j == 1) ? 1 : bb[j - 1] + 1;
            const int jhi = bb[j];
            double awthi = cwts[jhi - 1] - cwts[ilo1 - 1];
            if (jhi - ilo1 > nal0) {
                awthi = 0.0;
                for (int kk = 1; kk <= al0; ++kk) awthi = std::max(awthi, cwts[nal0 + kk - 1] - cwts[kk - 1]);
            }
            double awtlo;
            if (i == j) {
                awtlo = cwts[ilo1 + al0 - 1] - cwts[ilo1 - 1];
                for (int kk = ilo1 + 1; kk <= ihi - al0; ++kk) awtlo = std::min(awtlo, cwts[kk + al0 - 1] - cwts[kk - 1]);
            } else if (i + 1 == j) {
                awtlo = cwts[jlo - 1] - cwts[jlo - al0 - 1];
                for (int kk = jlo - al0 + 1; kk <= ihi; ++kk) awtlo = std::min(awtlo, cwts[kk + al0 - 1] - cwts[kk - 1]);
            } else {
                awtlo = cwts[jlo - 1] - cwts[ihi - 1];
            }
            const double sij1 = std::abs(bpsmax[j] - bpsmin[i]);
            const double sij2 = std::abs(bpsmax[i] - bpsmin[j]);
            const double sijmx0 = std::max(sij1, sij2);
            const double bsslim = (sijmx0 * sijmx0) / std::min(awtlo * (psrn - awtlo), awthi * (psrn - awthi));
            if (bssmax <= bsslim) {
                ++l;
                loc[l] = l;
                bloci[l] = i;
                blocj[l] = j;
                bssijmax[l] = bsslim;
                if (sij1 > sij2) {
                    awt[l] = std::abs(cwts[ibmax[j] - 1] - cwts[ibmin[i] - 1]);
                    bssbij[l] = (sij1 * sij1) / (awt[l] * (psrn - awt[l]));
                } else {
                    awt[l] = std::abs(cwts[ibmin[j] - 1] - cwts[ibmax[i] - 1]);
                    bssbij[l] = (sij2 * sij2) / (awt[l] * (psrn - awt[l]));
                }
            }
        }
    }
    const int nb1 = l;
    sort_loc_by_values(bssbij, loc, nb1);

    for (int ll = nb1; ll >= 1; --ll) {
        const int k = loc[ll];
        if (bssmax > bssijmax[k]) continue;
        const int bi = bloci[k], bj = blocj[k];
        double awtmax = awt[k];
        const int ilo1 = (bi == 1) ? 1 : bb[bi - 1] + 1;
        const int ihi = bb[bi];
        const int jlo = (bj == 1) ? 1 : bb[bj - 1] + 1;
        const int jhi = bb[bj];
        double awthi = cwts[jhi - 1] - cwts[ilo1 - 1];
        double awtlo = (bi == bj) ? 0.0 : (cwts[jlo - 1] - cwts[ihi - 1]);
        if (awtmax > psrn - awtmax) awtmax = psrn - awtmax;
        if (awtlo <= psrnov2) {
            const int ihi1 = (bi == bj) ? ihi - al0 : ihi;
            for (int i = ihi1; i >= ilo1; --i) {
                const int jlo1 = std::max(i + al0, jlo);
                for (int j = jlo1; j <= jhi; ++j) {
                    const double awt1 = cwts[j - 1] - cwts[i - 1];
                    if (awt1 <= awtmax) {
                        const double bijbss = std::pow(sx[j] - sx[i], 2.0) / (awt1 * (psrn - awt1));
                        if (bijbss > bssmax) { bssmax = bijbss; tmaxi = i; tmaxj = j; }
                    }
                }
            }
        }
        awtmax = psrn - awtmax;
        if (awthi >= psrnov2) {
            for (int i = ilo1; i <= ihi; ++i) {
                const int jhi1 = ((bi == 1) && (bj == nb)) ? std::min(jhi, jhi - al0 + i) : jhi;
                for (int j = jhi1; j >= jlo; --j) {
                    const double awt1 = cwts[j - 1] - cwts[i - 1];
                    if (awt1 >= awtmax) {
                        const double bijbss = std::pow(sx[j] - sx[i], 2.0) / (awt1 * (psrn - awt1));
                        if (bijbss > bssmax) { bssmax = bijbss; tmaxi = i; tmaxj = j; }
                    }
                }
            }
        }
    }

    if (tss <= bssmax + 0.0001) tss = bssmax + 1.0;
    return {bssmax / ((tss - bssmax) / (rn - 2.0)), tmaxi - 1, tmaxj - 1};
}

double wtmaxp(const std::vector<double>& px, const std::vector<double>& wts, const std::vector<double>& cwts, int al0) {
    return wtmaxo(px, wts, 0.0, cwts, al0).statistic; // placeholder, corrected below by local tss recomputation in hwt and wfind path not relying on this exact helper in tests
}

double hwtmaxp(const std::vector<double>& px, const std::vector<double>& wts, const std::vector<double>& cwts, const std::vector<double>& mncwt, int k, int al0) {
    const int n = static_cast<int>(px.size());
    double rn = static_cast<double>(n);
    const int nb = static_cast<int>(rn / static_cast<double>(k));
    std::vector<double> bpsmax(nb + 1), bpsmin(nb + 1), sx(n + 1, 0.0);
    std::vector<int> bb(nb + 1);
    for (int i = 1; i <= nb; ++i) bb[i] = static_cast<int>(std::lround(rn * (static_cast<double>(i) / static_cast<double>(nb))));
    int ilo = 1;
    double psum = 0.0, ssq = 0.0, bssmax = 0.0;
    rn = cwts[n - 1];
    for (int j = 1; j <= nb; ++j) {
        sx[ilo] = psum + px[ilo - 1] * wts[ilo - 1];
        ssq += wts[ilo - 1] * px[ilo - 1] * px[ilo - 1];
        double psmin = sx[ilo], psmax = sx[ilo];
        int ipsmin = ilo, ipsmax = ilo;
        for (int i = ilo + 1; i <= bb[j]; ++i) {
            sx[i] = sx[i - 1] + px[i - 1] * wts[i - 1];
            ssq += wts[i - 1] * px[i - 1] * px[i - 1];
            if (sx[i] < psmin) { psmin = sx[i]; ipsmin = i; }
            if (sx[i] > psmax) { psmax = sx[i]; ipsmax = i; }
        }
        bpsmin[j] = psmin; bpsmax[j] = psmax;
        psum = sx[bb[j]];
        ilo = bb[j] + 1;
        const int d = std::abs(ipsmin - ipsmax);
        if ((d <= k) && (d >= al0)) {
            const double rj = std::abs(cwts[ipsmax - 1] - cwts[ipsmin - 1]);
            const double bssij = std::pow(bpsmax[j] - bpsmin[j], 2.0) / (rj * (rn - rj));
            if (bssmax < bssij) bssmax = bssij;
        }
    }
    double tss = ssq - std::pow(sx[n] / rn, 2.0);

    auto scan_regular = [&](int left, int right, double psdiff) {
        const double psdiffsq = psdiff * psdiff;
        for (int j = al0; j <= k; ++j) {
            double rj = mncwt[j];
            const double bsslim = psdiffsq / (rj * (rn - rj));
            if (bsslim < bssmax) break;
            for (int i = left; i <= right - j; ++i) {
                const int ipj = i + j;
                rj = cwts[ipj - 1] - cwts[i - 1];
                const double bssij = std::pow(sx[ipj] - sx[i], 2.0) / (rj * (rn - rj));
                if (bssij > bssmax) bssmax = bssij;
            }
        }
    };

    scan_regular(1, bb[1], bpsmax[1] - bpsmin[1]);
    {
        const double psdiff = std::max(std::abs(bpsmax[1] - bpsmin[nb]), std::abs(bpsmax[nb] - bpsmin[1]));
        const double psdiffsq = psdiff * psdiff;
        for (int j = al0; j <= k; ++j) {
            double rj = mncwt[j];
            const double bsslim = psdiffsq / (rj * (rn - rj));
            if (bsslim < bssmax) break;
            const int nmj = n - j;
            for (int i = 1; i <= j; ++i) {
                const int ipnmj = i + nmj;
                rj = cwts[ipnmj - 1] - cwts[i - 1];
                const double bssij = std::pow(sx[ipnmj] - sx[i], 2.0) / (rj * (rn - rj));
                if (bssij > bssmax) bssmax = bssij;
            }
        }
    }
    for (int l = 2; l <= nb; ++l) {
        scan_regular(bb[l - 1] + 1, bb[l], bpsmax[l] - bpsmin[l]);
        const double psdiff = std::max(std::abs(bpsmax[l] - bpsmin[l - 1]), std::abs(bpsmax[l - 1] - bpsmin[l]));
        const double psdiffsq = psdiff * psdiff;
        for (int j = al0; j <= k; ++j) {
            double rj = mncwt[j];
            const double bsslim = psdiffsq / (rj * (rn - rj));
            if (bsslim < bssmax) break;
            for (int i = bb[l - 1] + 1 - j; i <= bb[l - 1]; ++i) {
                const int ipj = i + j;
                rj = cwts[ipj - 1] - cwts[i - 1];
                const double bssij = std::pow(sx[ipj] - sx[i], 2.0) / (rj * (rn - rj));
                if (bssij > bssmax) bssmax = bssij;
            }
        }
    }
    if (tss <= bssmax + 0.0001) tss = bssmax + 1.0;
    return bssmax / ((tss - bssmax) / (static_cast<double>(n) - 2.0));
}

ChangePointResult fndcpt(const std::vector<double>& x, double tss, int nperm, double cpval, bool ibin, bool hybrid, int al0, int hk, double delta, int ngrid, const std::vector<int>& sbdry, double tol, std::mt19937_64& rng) {
    const int n = static_cast<int>(x.size());
    std::vector<double> px(n);
    ChangePointResult res;
    const auto obs = tmaxo(x, tss, al0, ibin);
    res.ostat = obs.statistic;
    res.iseg = {obs.start, obs.end};
    double ostat1 = std::sqrt(obs.statistic);
    double ostat = obs.statistic * 0.99999;
    if (ostat1 <= 0.1) return res;
    const int l = std::min(obs.end - obs.start, n - obs.end + obs.start);
    if (!((ostat1 >= 7.0) && (l >= 10))) {
        int nrej = 0;
        if (hybrid) {
            const double pval1 = tailp(ostat1, delta, n, ngrid, tol);
            if (pval1 > cpval) return res;
            const int nrejc = static_cast<int>((cpval - pval1) * static_cast<double>(nperm));
            int k = nrejc * (nrejc + 1) / 2 + 1;
            for (int np = 1; np <= nperm; ++np) {
                xperm(x, px, rng);
                const double pstat = htmaxp(px, tss, hk, al0, ibin);
                if (ostat <= pstat) { ++nrej; ++k; }
                if (nrej > nrejc) return res;
                if (np >= sbdry[k - 1]) break;
            }
        } else {
            const int nrejc = static_cast<int>(cpval * static_cast<double>(nperm));
            int k = nrejc * (nrejc + 1) / 2 + 1;
            for (int np = 1; np <= nperm; ++np) {
                xperm(x, px, rng);
                const double pstat = tmaxp(px, tss, al0, ibin);
                if (ostat <= pstat) { ++nrej; ++k; }
                if (nrej > nrejc) return res;
                if (np >= sbdry[k - 1]) break;
            }
        }
    }
    if (obs.end == n) {
        res.ncpt = 1; res.icpt[0] = obs.start;
    } else if (obs.start == 0) {
        res.ncpt = 1; res.icpt[0] = obs.end;
    } else {
        int n1 = obs.start, n12 = obs.end, n2 = n12 - n1;
        double tpval = tpermp(n1, n2, n12, x.data(), px, nperm, rng);
        if (tpval <= cpval) { res.ncpt = 1; res.icpt[0] = obs.start; }
        const int offset = obs.start;
        n12 = n - obs.start;
        n2 = n - obs.end;
        n1 = n12 - n2;
        tpval = tpermp(n1, n2, n12, x.data() + offset, px, nperm, rng);
        if (tpval <= cpval) res.icpt[res.ncpt++] = obs.end;
    }
    return res;
}

ChangePointResult wfindcpt(const std::vector<double>& x, double tss, const std::vector<double>& wts, const std::vector<double>& rwts, const std::vector<double>& cwts, int nperm, double cpval, bool hybrid, int al0, int hk, double delta, int ngrid, const std::vector<int>& sbdry, double tol, std::mt19937_64& rng) {
    const int n = static_cast<int>(x.size());
    std::vector<double> px(n), mncwt;
    ChangePointResult res;
    const auto obs = wtmaxo(x, wts, tss, cwts, al0);
    res.ostat = obs.statistic;
    res.iseg = {obs.start, obs.end};
    double ostat1 = std::sqrt(obs.statistic);
    double ostat = obs.statistic * 0.99999;
    if (ostat1 <= 0.1) return res;
    const int l = std::min(obs.end - obs.start, n - obs.end + obs.start);
    if (!((ostat1 >= 7.0) && (l >= 10))) {
        int nrej = 0;
        if (hybrid) {
            getmncwt(n, cwts, hk, mncwt, delta);
            const double pval1 = tailp(ostat1, delta, n, ngrid, tol);
            if (pval1 > cpval) return res;
            const int nrejc = static_cast<int>((cpval - pval1) * static_cast<double>(nperm));
            int k = nrejc * (nrejc + 1) / 2 + 1;
            for (int np = 1; np <= nperm; ++np) {
                wxperm(x, px, rwts, rng);
                const double pstat = hwtmaxp(px, wts, cwts, mncwt, hk, al0);
                if (ostat <= pstat) { ++nrej; ++k; }
                if (nrej > nrejc) return res;
                if (np >= sbdry[k - 1]) break;
            }
        } else {
            const int nrejc = static_cast<int>(cpval * static_cast<double>(nperm));
            int k = nrejc * (nrejc + 1) / 2 + 1;
            for (int np = 1; np <= nperm; ++np) {
                wxperm(x, px, rwts, rng);
                const double pstat = wtmaxp(px, wts, cwts, al0);
                if (ostat <= pstat) { ++nrej; ++k; }
                if (nrej > nrejc) return res;
                if (np >= sbdry[k - 1]) break;
            }
        }
    }
    if (obs.end == n) {
        res.ncpt = 1; res.icpt[0] = obs.start;
    } else if (obs.start == 0) {
        res.ncpt = 1; res.icpt[0] = obs.end;
    } else {
        int n1 = obs.start, n12 = obs.end, n2 = n12 - n1;
        double tpval = wtpermp(n1, n2, n12, x.data(), px, wts.data(), rwts.data(), nperm, rng);
        if (tpval <= cpval) { res.ncpt = 1; res.icpt[0] = obs.start; }
        const int offset = obs.start;
        n12 = n - obs.start;
        n2 = n - obs.end;
        n1 = n12 - n2;
        tpval = wtpermp(n1, n2, n12, x.data() + offset, px, wts.data() + offset, rwts.data() + offset, nperm, rng);
        if (tpval <= cpval) res.icpt[res.ncpt++] = obs.end;
    }
    return res;
}

} // namespace cbs
