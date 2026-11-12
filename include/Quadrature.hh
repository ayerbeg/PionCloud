#pragma once
// ============================================================
//  Quadrature.hh   (new in v1.4)
//  Gauss-Legendre nodes/weights and the exact |k|/angle region of the
//  kbin scan in the variable u = kT^2.
//
//  Used by Calculator::kbinF2Exact() (kbin_method = exact).  See the
//  comments there and CHANGELOG.md ("kbin: exact integration region").
// ============================================================

#include "PhysicsParams.hh"
#include <cmath>
#include <utility>
#include <vector>

namespace PionCloud {

// ----------------------------------------------------------
//  Gauss-Legendre rule on [-1, 1] with n nodes (Newton iteration on
//  P_n, as in Numerical Recipes "gauleg").  Exact for polynomials of
//  degree <= 2n-1.
// ----------------------------------------------------------
struct GaussLegendre {
    std::vector<double> x, w;
    explicit GaussLegendre(int n = 16) : x(n), w(n) {
        const int m = (n + 1) / 2;
        for (int i = 0; i < m; ++i) {
            double z = std::cos(PI * (i + 0.75) / (n + 0.5)), z1, pp;
            do {
                double p1 = 1.0, p2 = 0.0;
                for (int j = 1; j <= n; ++j) {
                    const double p3 = p2; p2 = p1;
                    p1 = ((2.0*j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
                }
                pp = n * (z * p1 - p2) / (z*z - 1.0);
                z1 = z; z = z1 - p1 / pp;
            } while (std::fabs(z - z1) > 1e-15);
            x[i] = -z;  x[n-1-i] = z;
            w[i] = w[n-1-i] = 2.0 / ((1.0 - z*z) * pp*pp);
        }
    }
    // integral of f over [a, b]
    template <class F> double integrate(F &&f, double a, double b) const {
        const double h = 0.5 * (b - a), c = 0.5 * (b + a);
        double s = 0.0;
        for (size_t i = 0; i < x.size(); ++i) s += w[i] * f(c + h * x[i]);
        return s * h;
    }
};

// ----------------------------------------------------------
//  Allowed region in u = kT^2 at fixed meson fraction y.
//
//  With A = 2 mN (1-y) and c = (1-(1-y)^2) mN^2 = y(2-y) mN^2:
//     |k|^2(u)   = u + (u + c)^2 / A^2          (monotonic in u)
//     cosphi(u)  = (u + c) / (A |k|)            (> 0 always)
//  (the kMag/cosPhi formulas of f2_pi_sub.f90 written in u).
//
//  Conditions:
//   kmin <= |k| <= kmax :  u between the roots of
//                          u^2 + (2c + A^2) u + c^2 - A^2 K^2 = 0
//   cosphi <= Cmax      :  (1-Cmax^2)(u+c)^2 - Cmax^2 A^2 u <= 0
//                          -> u inside [r1, r2]  (empty if no real roots)
//   cosphi >= Cmin      :  (1-Cmin^2)(u+c)^2 - Cmin^2 A^2 u >= 0
//                          -> u outside (s1, s2)
//   kT <= kTmax         :  u <= kTmax^2  (same 10 GeV as the grid method)
//  Result: at most two intervals.
// ----------------------------------------------------------
inline std::vector<std::pair<double,double>>
kbinRegionU(double y, double kmin, double kmax,
            double cmin, double cmax, double kTmax = 10.0)
{
    std::vector<std::pair<double,double>> out;
    if (y <= 0.0 || y >= 1.0) return out;
    const double A  = 2.0 * mN * (1.0 - y), A2 = A*A;
    const double c  = y * (2.0 - y) * mN*mN;
    const double b  = 2.0*c + A2;
    // u at which |k| = K (largest root; <= 0 if |k|(u=0) >= K)
    auto uOfK = [&](double K) {
        const double cc = c*c - A2*K*K;            // value at u = 0
        if (cc >= 0.0) return 0.0;                 // |k|(0) >= K
        return -2.0*cc / (b + std::sqrt(b*b - 4.0*cc));
    };
    if (c / A > kmax) return out;                  // |k|(u=0) = c/A already > kmax
    double lo = uOfK(kmin), hi = std::min(uOfK(kmax), kTmax*kTmax);
    if (hi <= lo) return out;

    auto roots = [&](double C, double &r1, double &r2) -> bool {
        const double a2 = 1.0 - C*C;
        const double b2 = 2.0*a2*c - C*C*A2;
        const double c2 = a2*c*c;
        const double d  = b2*b2 - 4.0*a2*c2;
        if (d <= 0.0 || a2 <= 0.0) return false;
        const double sq = std::sqrt(d);
        r1 = (-b2 - sq) / (2.0*a2);
        r2 = (-b2 + sq) / (2.0*a2);
        return true;
    };
    // cosphi <= cmax
    if (cmax < 1.0) {
        double r1, r2;
        if (!roots(cmax, r1, r2)) return out;      // cosphi > cmax everywhere
        lo = std::max(lo, r1); hi = std::min(hi, r2);
        if (hi <= lo) return out;
    }
    // cosphi >= cmin : remove (s1, s2)
    double s1, s2;
    if (cmin > 0.0 && roots(cmin, s1, s2)) {
        if (s1 > lo) out.emplace_back(lo, std::min(hi, s1));
        if (s2 < hi) out.emplace_back(std::max(lo, s2), hi);
        std::vector<std::pair<double,double>> ok;
        for (auto &iv : out) if (iv.second > iv.first) ok.push_back(iv);
        return ok;
    }
    out.emplace_back(lo, hi);
    return out;
}

// Largest meson fraction y for which |k| <= kmax is reachable
// (|k| at kT = 0 equals kmax):  mN y^2 - 2(mN + kmax) y + 2 kmax = 0.
inline double kbinYmax(double kmax)
{
    const double a = mN + kmax;
    return (a - std::sqrt(a*a - 2.0*mN*kmax)) / mN;
}

} // namespace PionCloud
