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
//  Convenience: F2pi from GRV92
//  F2pi = (5/9) * (xVpi + 2*xSpi)   [quark model charge weighting]
// ----------------------------------------------------------
inline double F2piGRV(double x, double Q2)
{
    double xV, xS;
    GRV92(x, Q2, xV, xS);
    return (5.0/9.0) * (xV + 2.0*xS);
}

// ----------------------------------------------------------
//  GRV99  -- LO pion valence + sea distributions (1999 update)
//
//  Translated from SUBROUTINE GRV99(x,Q2,xVpi,xSpi) in
//  f2_pi_sub.f90 lines 766-823.
//
//  Despite the Fortran comment citing Z.Phys.C53(1992) this
//  is the 1999 parametrisation: Q02=0.26 GeV^2, L=0.204 GeV,
//  and different coefficients from GRV92 (Q02=0.25, L=0.232).
//
//  Valid for 0.25 < Q2 < 1e8 GeV^2, 1e-5 < x < 1.
//
//  NOTE on the Fortran source: the variable BB is used but not
//  explicitly declared in the REAL*8 list (implicit typing).
//  It corresponds to the b coefficient in the valence term.
//  Reproduced here as the local variable 'BB'.
// ----------------------------------------------------------
inline void GRV99(double x, double Q2, double &xVpi, double &xSpi)
{
    xVpi = 0.0;
    xSpi = 0.0;

    if (x  < 1.0e-5) x  = 1.01e-5;
    if (Q2 < 0.25)   Q2 = 0.2501;

    const double Q02 = 0.26;    // mu^2_LO
    const double L   = 0.204;
    if (Q2 <= Q02) return;

    const double s = std::log( std::log(Q2/(L*L)) / std::log(Q02/(L*L)) );

    // ---- valence distribution ----
    const double Nv = 1.212 + 0.498*s - 0.009*s*s;
    const double a  = 0.517 - 0.020*s;
    const double AA = -0.037 - 0.578*s;
    const double BB =  0.241 + 0.251*s;   // 'BB' from Fortran (implicitly typed)
    const double D  =  0.383 + 0.624*s;

    if (x < 1.0)
        xVpi = Nv * std::pow(x, a)
             * (1.0 + AA*std::sqrt(x) + BB*x)
             * std::pow(1.0 - x, D);

    // ---- sea distribution (SU(3) symmetric) ----
    const double alpha = 0.823;
    const double as_   = 1.036 - 0.709*s;
    const double AAs   = -1.245 + 0.713*s;
    const double Bs    =  5.580 - 1.281*s;
    const double Ds    =  2.746 - 0.191*s;
    const double E_    =  5.101 + 1.294*s;
    const double Epr   =  4.854 - 0.437*s;
    const double beta  =  0.650;

    if (x < 1.0 && s > 0.0)
        xSpi = std::pow(s, alpha) / std::pow(std::log(1.0/x), as_)
             * (1.0 + AAs*std::sqrt(x) + Bs*x)
             * std::pow(1.0 - x, Ds)
             * std::exp(-E_ + std::sqrt(Epr * std::pow(s, beta) * std::log(1.0/x)));
}

// ----------------------------------------------------------
//  Convenience: F2pi from GRV99
// ----------------------------------------------------------
inline double F2piGRV99(double x, double Q2)
{
    double xV, xS;
    GRV99(x, Q2, xV, xS);
    return (5.0/9.0) * (xV + 2.0*xS);
}

} // namespace PionCloud
