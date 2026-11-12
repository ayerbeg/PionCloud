#pragma once
// ============================================================
//  GRV.hh
//  LO pion parton distributions
//  Gluck, Reya, Vogt: Z.Phys. C53 (1992) 651, appendix 1
//
//  Translated from SUBROUTINE GRV(x,Q2,xVpi,xSpi)
//  Valid for 0.25 < Q2 < 1e8 GeV^2, 1e-5 < x < 1
// ============================================================

#include <cmath>
#include <algorithm>

namespace PionCloud {

// ----------------------------------------------------------
//  GRV92  -- LO pion valence + sea distributions
//  Returns xVpi (valence) and xSpi (sea).
// ----------------------------------------------------------
inline void GRV92(double x, double Q2, double &xVpi, double &xSpi)
{
    xVpi = 0.0;
    xSpi = 0.0;

    if (x   < 1.0e-5) x  = 1.01e-5;
    if (Q2  < 0.25)   Q2 = 0.2501;

    const double Q02 = 0.25;
    const double L   = 0.232;
    if (Q2 <= Q02) return;

    double s = std::log( std::log(Q2/(L*L)) / std::log(Q02/(L*L)) );

    // --- Valence distribution ---
    double Nv = 0.519 + 0.180*s - 0.011*s*s;
    double a  = 0.499 - 0.027*s;
    double AA = 0.381 - 0.419*s;
    double D  = 0.367 + 0.563*s;
    if (x < 1.0)
        xVpi = Nv * std::pow(x, a) * (1.0 + AA*std::sqrt(x)) * std::pow(1.0-x, D);

    // --- Sea distribution (SU(3) symmetric) ---
    double alpha = 0.55;
    double as_   = 2.538 - 0.763*s;
    double AAs   = -0.748;
    double Bs    = 0.313 + 0.935*s;
    double Ds    = 3.359;
    double E_    = 4.433 + 1.301*s;
    double Epr   = 9.30  - 0.887*s;
    double beta  = 0.56;
    if (x < 1.0)
        xSpi = std::pow(s, alpha) / std::pow(std::log(1.0/x), as_)
             * (1.0 + AAs*std::sqrt(x) + Bs*x)
             * std::pow(1.0-x, Ds)
             * std::exp(-E_ + std::sqrt(Epr * std::pow(s, beta) * std::log(1.0/x)));
}

// ----------------------------------------------------------
//  Convenience: F2pi from GRV92 (LO, quark-parton model)
//
//  For pi+ = (u dbar) with per-flavour valence V = u_v = dbar_v and
//  an SU(3)-symmetric sea S = ubar = u_sea = d = dbar_sea = s = sbar:
//
//    F2pi = x [ 4/9 (u+ubar) + 1/9 (d+dbar) + 1/9 (s+sbar) ]
//         = 5/9 (xV + 2 xS)            <- u,d quarks   (Fortran, default)
//         + 2/9 xS                     <- strange quarks (optional)
//
//  The original Fortran (3Var_*.f, f2_pi_sub.f90) keeps only the
//  u,d term; include_strange=true adds the strange-sea term.
//  The same expression holds for pi0 and pi- at LO (isospin).
// ----------------------------------------------------------
inline double F2piGRV(double x, double Q2, bool include_strange = false)
{
    double xV, xS;
    GRV92(x, Q2, xV, xS);
    double f2 = (5.0/9.0) * (xV + 2.0*xS);
    if (include_strange) f2 += (2.0/9.0) * xS;   // SU(3) sea: sbar = S
    return f2;
}

// ============================================================
//  GRS99 -- Gluck, Reya, Schienbein, "Pionic parton distributions
//           revisited", Eur. Phys. J. C10 (1999) 313,
//           arXiv:hep-ph/9903288, Appendix A (LO).
//
//  (Called "GRV99" in this code and in f2_pi_sub.f90 for historical
//   reasons; the parametrisation is GRS99.  Its header comment in the
//   Fortran wrongly cites Z.Phys.C53 (1992) 651.)
//
//  CORRECTED in v1.4 (see CHANGELOG.md, "GRV99 fix").  The v1.3
//  translation of the Fortran routine GRV99 (f2_pi_sub.f90:766-823,
//  whose call was commented out and never used there) had three errors:
//
//   1. VALENCE NORMALISATION.  GRS99 define (text after their Eq. 2)
//        v^pi = u_v^pi + dbar_v^pi        (BOTH valence quarks)
//      so  int v^pi dx = 2.  GRV92's xV is ONE valence flavour
//      (int V dx = 1).  v1.3 passed GRS99's total v^pi into the per-
//      flavour formula 5/9(xV+2xS) -> F2pi^val twice too large
//      (F2^(piN) +107% on the TDIS Fig. 37 configuration).
//      Now: xVpi = v^pi / 2.
//   2. WRONG SEA.  v1.3 used the coefficients of Eq. (A.7)-(A.8),
//      which are the STRANGE sea s^pi = sbar^pi, as the light sea.
//      Now: light sea qbar^pi = (ubar+dbar)/2 from Eq. (A.4)+(A.6);
//      strange sea returned separately.
//   3. SIGN.  N = 1.212 + 0.498 s + 0.009 s^2 (Eq. A.3); v1.3 had
//      -0.009 s^2.
//
//  Checks (tests/test_all.cc, test_GRS99_sum_rules):
//    int v^pi dx = 2.00 and the momentum sum
//    int x (v + 4 qbar + 2 sbar + g) dx = 1.00 +- 0.005
//    for Q2 = 1..100 GeV^2  -- both hold only with the corrections.
//
//  Validity (paper): 0.5 <~ Q2 <~ 1e5 GeV^2, 1e-5 < x < 1.
//  Below mu^2_LO = 0.26 GeV^2 the input scale (s = 0) is used.
// ============================================================
struct GRS99Partons {
    double xv    = 0.0;   // x v^pi    = x(u_v + dbar_v), TOTAL valence (int v = 2)
    double xqbar = 0.0;   // x qbar^pi = x(ubar + dbar)/2, light sea per flavour
    double xsbar = 0.0;   // x sbar^pi = x s^pi, strange sea per flavour
    double xg    = 0.0;   // x g^pi
};

namespace detail {
// Common functional form of Eq. (A.4) for gluon and light sea:
//   x w = [ x^a (A + B sqrt(x) + C x) (ln 1/x)^b
//           + s^alpha exp(-E + sqrt(E' s^beta ln 1/x)) ] (1-x)^D
inline double grs99_w(double x, double s, double alpha, double beta,
                      double a, double b, double A, double B, double C,
                      double D, double E, double Ep)
{
    const double L = std::log(1.0/x);
    return ( std::pow(x, a) * (A + B*std::sqrt(x) + C*x) * std::pow(L, b)
           + std::pow(s, alpha) * std::exp(-E + std::sqrt(Ep * std::pow(s, beta) * L)) )
         * std::pow(1.0 - x, D);
}
} // namespace detail

inline GRS99Partons GRS99(double x, double Q2)
{
    GRS99Partons f;
    if (x >= 1.0) return f;
    if (x < 1.0e-5) x = 1.01e-5;

    const double mu2 = 0.26;     // mu^2_LO  [GeV^2]
    const double L   = 0.204;    // Lambda_LO [GeV], Eq. (A.1)
    if (Q2 < mu2) Q2 = mu2;      // freeze at the input scale (s = 0)
    const double s  = std::log( std::log(Q2/(L*L)) / std::log(mu2/(L*L)) );
    const double rs = std::sqrt(s);

    // ---- valence, Eqs. (A.2)-(A.3) ----
    {
        const double N = 1.212 + 0.498*s + 0.009*s*s;
        const double a = 0.517 - 0.020*s;
        const double A = -0.037 - 0.578*s;
        const double B =  0.241 + 0.251*s;
        const double D =  0.383 + 0.624*s;
        f.xv = N * std::pow(x, a) * (1.0 + A*std::sqrt(x) + B*x) * std::pow(1.0 - x, D);
    }
    // ---- gluon, Eqs. (A.4)-(A.5) ----
    f.xg = detail::grs99_w(x, s, 0.504, 0.226,
                           2.251 - 1.339*rs, 0.0,
                           2.668 - 1.265*s + 0.156*s*s, -1.839 + 0.386*s,
                           -1.014 + 0.920*s - 0.101*s*s,
                           -0.077 + 1.466*s, 1.245 + 1.833*s, 0.510 + 3.844*s);
    // ---- light sea, Eqs. (A.4) and (A.6) ----
    f.xqbar = detail::grs99_w(x, s, 1.147, 1.241,
                              0.309 - 0.134*rs, 0.893 - 0.264*rs,
                              0.219 - 0.054*s, -0.593 + 0.240*s, 1.100 - 0.452*s,
                              3.526 + 0.491*s, 4.521 + 1.583*s, 3.102);
    // ---- strange sea, Eqs. (A.7)-(A.8)  (zero at the input scale, s = 0) ----
    if (s > 0.0) {
        const double Ls = std::log(1.0/x);
        const double a  = 1.036 - 0.709*s;
        f.xsbar = std::pow(s, 0.823) / std::pow(Ls, a)
                * (1.0 + (-1.245 + 0.713*s)*std::sqrt(x) + (5.580 - 1.281*s)*x)
                * std::pow(1.0 - x, 2.746 - 0.191*s)
                * std::exp(-(5.101 + 1.294*s) + std::sqrt((4.854 - 0.437*s)*std::pow(s, 0.650)*Ls));
    }
    return f;
}

// ----------------------------------------------------------
//  GRV99 -- per-flavour interface, SAME CONVENTION AS GRV92:
//    xVpi = valence of one flavour  (= x v^pi / 2)
//    xSpi = light sea of one flavour (= x qbar^pi)
//  so that F2pi = 5/9 (xVpi + 2 xSpi) holds for both sets.
// ----------------------------------------------------------
inline void GRV99(double x, double Q2, double &xVpi, double &xSpi)
{
    const GRS99Partons f = GRS99(x, Q2);
    xVpi = 0.5 * f.xv;
    xSpi = f.xqbar;
}

// ----------------------------------------------------------
//  Convenience: F2pi from GRS99 (LO).  include_strange adds
//  2/9 x sbar (GRS99 has an SU(3)-broken strange sea).
// ----------------------------------------------------------
inline double F2piGRV99(double x, double Q2, bool include_strange = false)
{
    const GRS99Partons f = GRS99(x, Q2);
    double f2 = (5.0/9.0) * (0.5*f.xv + 2.0*f.xqbar);
    if (include_strange) f2 += (2.0/9.0) * f.xsbar;
    return f2;
}

} // namespace PionCloud
