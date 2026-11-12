// ============================================================
//  Calculator.cc
//  Pion Cloud Model -- integration engine implementation
//
//  All loops use integer indices to avoid floating-point
//  accumulation errors in the loop variable.  The physical
//  coordinate is computed as:
//      coord = min + i * step
//  which is exact for each step and avoids the classic
//  "one extra iteration" bug.
// ============================================================

#include "Calculator.hh"
#include "SplittingFunctions.hh"
#include "GRV.hh"
#include "Kinematics.hh"

#ifdef HAVE_CTEQ6
#include "CTEQ6.hh"
#endif

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#ifndef _WIN32
#  include <unistd.h>   // chdir, getcwd
#else
#  include <direct.h>
#endif

namespace PionCloud {

Calculator::Calculator(const RunParams &p) : p_(p)
{
    // ---- apply precision preset (overrides ny/nkT if not CUSTOM) ----
    switch (p_.precision) {
    case PRECISION_FAST:
        p_.ny  = 20;   p_.nkT = 500;   break;
    case PRECISION_NORMAL:
        p_.ny  = 50;   p_.nkT = 2000;  break;
    case PRECISION_PUBLICATION:
        p_.ny  = 100;  p_.nkT = 50000; break;
    case PRECISION_CUSTOM:
    default:
        break;  // use ny/nkT as set by user
    }

    // ---- channel / charge-state consistency (v1.4) ----
    // The Delta channels have zero coupling for dis = neutral in the
    // Fortran (fypiD, f_RhoDel); v1.3 then returned F2 = 0 silently.
    if ((p_.flag == FLAG_PI_DELTA || p_.flag == FLAG_RHO_DELTA) &&
        p_.dis == DIS_NEUTRAL)
        throw std::invalid_argument(
            "Calculator: pi-Delta / rho-Delta are only implemented for "
            "dissociation = charge_exchange (sum of charge states); "
            "dissociation = neutral would give F2 = 0.");

    // ---- regulator / cut-off plausibility (v1.4, warning only) ----
    // Published cut-offs differ by a factor ~5 between regulators
    // (McKenney et al. PRD 93 054011 Table I: s-exp 1.33, t-exp 0.85,
    // Pauli-Villars 0.27 GeV; HSS s-exp 1.56 GeV).  A cut-off far from the
    // range of the chosen regulator usually means a forgotten form_factor
    // or lambda setting (as in the v1.3 default PV + 1.33 GeV).
    {
        const bool   isDelta = (p_.flag == FLAG_PI_DELTA || p_.flag == FLAG_RHO_DELTA);
        const double lam     = isDelta ? p_.Ld : p_.L;   // the cut-off this channel uses
        auto warn = [&](const char *ff, double lo, double hi) {
                if (lam < lo || lam > hi)
                    std::cerr << "[Calculator] WARNING: form factor " << ff
                              << " with Lambda = " << lam << " GeV is outside the"
                              << " published range [" << lo << ", " << hi
                              << "] GeV for this regulator -- check form_factor/lambda.\n";
        };
        if (p_.typ == FF_PAULI_VILLAR) warn("pauli_villar", 0.15, 0.6);
        if (p_.typ == FF_EXP_S)        warn("exp_s",        0.9,  2.0);
        if (p_.typ == FF_EXP_T)        warn("exp_t",        0.5,  1.5);
    }

    // ---- inclusive mode: remove all kinematic cuts ----
    if (p_.inclusive_mode) {
        p_.cosph_min = 0.0;
        p_.cosph_max = 1.0;
        p_.kmin      = 0.001;
        p_.kmax      = 10.0;
    }

#ifdef HAVE_CTEQ6
    // If the user specified a table directory, chdir() there before
    // calling CTEQ6::init() (SetCtq6 opens files by name in the cwd),
    // then restore the working directory afterwards.
    {
        char savedCwd[4096] = {0};
#ifdef _WIN32
        _getcwd(savedCwd, sizeof(savedCwd));
#else
        if (getcwd(savedCwd, sizeof(savedCwd)) == nullptr) savedCwd[0] = '\0';
#endif
        if (!p_.cteq6_table_dir.empty()) {
            if (chdir(p_.cteq6_table_dir.c_str()) != 0)
                std::cerr << "[Calculator] WARNING: could not chdir to '"
                          << p_.cteq6_table_dir << "' -- CTEQ6 may not find tables.\n";
        }
        // v1.4: make sure the table is there BEFORE calling the Fortran.
        // If SetCtq6 cannot open it, it executes a Fortran STOP, which
        // terminates the whole host process with exit status 0 -- a silent
        // "success".  This hid the failure of most unit tests in v1.3
        // (ctest reported Passed) and would kill any program embedding
        // the library.  Throwing here lets the caller handle it.
        {
            std::ifstream tbl("cteq6m.tbl");
            if (!tbl.good()) {
                std::string where = p_.cteq6_table_dir.empty()
                                  ? std::string("the current working directory")
                                  : "'" + p_.cteq6_table_dir + "'";
                if (!p_.cteq6_table_dir.empty() && savedCwd[0] != '\0') chdir(savedCwd);
                throw std::runtime_error(
                    "Calculator: CTEQ6 table 'cteq6m.tbl' not found in " + where +
                    ".  Run from the build directory, copy the CTEQ6/*.tbl files, or "
                    "set cteq6_table_dir.");
            }
        }
        CTEQ6::init(CTEQ6::CTEQ6M);
        if (!p_.cteq6_table_dir.empty() && savedCwd[0] != '\0')
            chdir(savedCwd);   // restore
        if (p_.verbose)
            std::cout << "[Calculator] CTEQ6 initialised (CTEQ6M, iset=1, NLO MSbar).\n";
    }
#else
    if (p_.verbose)
        std::cout << "[Calculator] Built without CTEQ6 -- F2n and ratio will be zero.\n";
#endif
}

// ----------------------------------------------------------
//  integrand -- evaluates the full integrand at one phase-
//               space point (x, Q2, theta_h, k).
//
//  Steps:
//   1. (k, theta_h) -> (y, kT) coordinate transform
//   2. Kinematic guard (y out of range, x/y >= 1)
//   3. Phase-space Jacobian rho_PS
//   4. Meson-baryon splitting function
//   5. Pion structure function F2pi from GRV LO PDFs
// ----------------------------------------------------------
double Calculator::integrand(double x, double Q2,
                             double theta_h, double k) const
{
    double y, kT;
    kThetaToYkT(k, theta_h, y, kT);

    if (y <= 0.0 || y >= 1.0)  return 0.0;
    if (x <= 0.0 || x / y >= 1.0) return 0.0;
    if (Q2 <= 0.0)             return 0.0;

    double rho_PS   = phaseSpace(k, theta_h);
    double f2pi_grv = f2pi(x / y, Q2);
    double split    = splitting(y, kT);

    return rho_PS * split * f2pi_grv;
}

// ----------------------------------------------------------
//  splitting -- f_MB(y, kT) of the selected channel.
//  The ONLY place where the channel is dispatched; to add a new
//  meson-baryon channel, add its function to SplittingFunctions.hh,
//  a ContributionFlag value, and one case here.
// ----------------------------------------------------------
double Calculator::splitting(double y, double kT) const
{
    switch (p_.flag) {
    case FLAG_PION:      return fypiN   (y, kT, p_.L,  p_.typ, p_.dis);
    case FLAG_RHO:       return f_rhoN  (y, kT, p_.L,  p_.typ, p_.dis);
    case FLAG_PI_DELTA:  return fypiD   (y, kT, p_.Ld, p_.typ, p_.dis);
    case FLAG_RHO_DELTA: return f_RhoDel(y, kT, p_.Ld, p_.typ, p_.dis);
    default:
        throw std::invalid_argument("Calculator: Unknown ContributionFlag");
    }
}

// ----------------------------------------------------------
//  f2pi -- pion structure function at (x_pi, Q2).
//  The ONLY place where the pion PDF set is dispatched.
// ----------------------------------------------------------
double Calculator::f2pi(double xpi, double Q2) const
{
    return (p_.grv_version == GRV_99)
         ? PionCloud::F2piGRV99(xpi, Q2, p_.pion_strange)
         : PionCloud::F2piGRV  (xpi, Q2, p_.pion_strange);
}

// ----------------------------------------------------------
//  runXScan -- mirrors 3Var_x.f
//
//  Output loop:  x          [xmin, xmax]          (nx+1 points)
//  Middle loop:  theta_h    [thmin_deg, thmax_deg] (nth_int+1 pts, Simpson)
//  Inner loop:   k          [kmin, kmax]           (nk+1 pts, Simpson)
//
//  The outer x loop produces the output array.  theta_h and k
//  are both numerically integrated at each x with Simpson's
//  composite rule.
// ----------------------------------------------------------
void Calculator::runXScan(ScanResults &out) const
{
    out.xScan.clear();
    out.xScan.reserve(p_.nx + 1);

    const double xstep  = (p_.xmax - p_.xmin) / p_.nx;

    const double thmin  = (p_.thmin_deg / 180.0) * PI;
    const double thmax  = (p_.thmax_deg / 180.0) * PI;
    const double thstep = (thmax - thmin) / p_.nth_int;

    const double kstep  = (p_.kmax - p_.kmin) / p_.nk;

    double F2pi0 = 0.0;   // first moment ∫ F2piK(x) dx

    std::cout << std::scientific << std::setprecision(5);

    for (int ix = 0; ix <= p_.nx; ++ix) {
        const double x  = p_.xmin + ix * xstep;
        const double Q2 = calcQ2(p_.E, x, p_.theta_e);

        // ---- theta_h integration (Simpson composite) ----
        double F2piK = 0.0;
        for (int ith = 0; ith <= p_.nth_int; ++ith) {
            const double th  = thmin + ith * thstep;
            const double wth = simpsonWeight(ith, p_.nth_int);

            // ---- |k| integration (Simpson composite) ----
            double F2pi_k = 0.0;
            for (int ik = 0; ik <= p_.nk; ++ik) {
                const double k  = p_.kmin + ik * kstep;
                const double wk = simpsonWeight(ik, p_.nk);
                F2pi_k += wk * integrand(x, Q2, th, k);
            }
            F2pi_k *= (kstep / 3.0);

            F2piK += wth * F2pi_k;
        }
        F2piK *= (thstep / 3.0);

        F2pi0 += simpsonWeight(ix, p_.nx) * F2piK;

        XScanPoint pt;
        pt.x           = x;
        pt.Q2          = Q2;
        pt.F2piK       = F2piK;
        pt.F2piMoment  = 0.0;
        pt.F2n         = 0.0;
        pt.ratio       = 0.0;

#ifdef HAVE_CTEQ6
        pt.F2n = (p_.nucleon == NUCLEON_PROTON)
               ? CTEQ6::F2proton (x, Q2, false)
               : CTEQ6::F2neutron(x, Q2, false);
        if (pt.F2n > 0.0)
            pt.ratio = pt.F2piK / pt.F2n;
#endif
        out.xScan.push_back(pt);

        if (p_.verbose)
            std::cout << "  [x-scan]  x=" << x
                      << "  Q2=" << Q2
                      << "  F2piK=" << F2piK
#ifdef HAVE_CTEQ6
                      << "  F2n=" << pt.F2n
                      << "  ratio=" << pt.ratio
#endif
                      << "\n";
    }

    out.F2pi0_x = (xstep / 3.0) * F2pi0;

    // Running moment using composite Simpson (consistent with first moment)
    {
        int np = static_cast<int>(out.xScan.size());
        double cumulative = 0.0;
        for (int i = 0; i < np; ++i) {
            cumulative += simpsonWeight(i, np - 1) * out.xScan[i].F2piK;
            out.xScan[i].F2piMoment = (xstep / 3.0) * cumulative;
        }
    }

    if (p_.verbose)
        std::cout << "  [x-scan]  First moment of F2piK(x) = "
                  << out.F2pi0_x << "\n";
}

// ----------------------------------------------------------
//  runThetaScan -- mirrors 3Var_theta.f
//
//  Output loop:  theta_h    [thmin_deg, thmax_deg]  (nth+1 points)
//  Middle loop:  x          [xint_th_start, xmax_th_fixed] (nx_int+1, Simpson)
//  Inner loop:   k          [kmin_th_fixed, kmax_th_fixed] (nk+1,     Simpson)
//
//  The x and k integration ranges cover the full accessible
//  phase space, as in 3Var_theta.f.  They are different from
//  the narrow BoNuS window used in runXScan; see Calculator.hh
//  for the fixed-range constants.
// ----------------------------------------------------------
void Calculator::runThetaScan(ScanResults &out) const
{
    out.thetaScan.clear();
    out.thetaScan.reserve(p_.nth + 1);

    const double thmin  = (p_.thmin_deg / 180.0) * PI;
    const double thmax  = (p_.thmax_deg / 180.0) * PI;
    const double thstep = (thmax - thmin) / p_.nth;

    // x and k integration ranges from RunParams
    // Left panel defaults: xmin_th~0, xmax_th=0.6, kmin_th~0, kmax_th=0.5
    // Right panel:         xmin_th=0.05,            kmin_th=0.060, kmax_th=0.250
    const double xmin_th  = p_.xmin_th;
    const double xmax_th  = p_.xmax_th;
    const double xstep_th = (xmax_th - xmin_th) / p_.nx_int;

    const double kmin_th  = p_.kmin_th;
    const double kmax_th  = p_.kmax_th;
    const double kstep_th = (kmax_th - kmin_th) / p_.nk;

    double F2pi0 = 0.0;   // first moment ∫ F2piK(θ) dθ

    std::cout << std::scientific << std::setprecision(5);

    for (int ith = 0; ith <= p_.nth; ++ith) {
        const double th  = thmin + ith * thstep;
        double F2piK = 0.0;

        // ---- x integration (Simpson composite) ----
        for (int ix = 0; ix <= p_.nx_int; ++ix) {
            const double x  = xmin_th + ix * xstep_th;
            const double wx = simpsonWeight(ix, p_.nx_int);
            const double Q2 = calcQ2(p_.E, x, p_.theta_e);

            // ---- |k| integration (Simpson composite) ----
            double F2pi_k = 0.0;
            for (int ik = 0; ik <= p_.nk; ++ik) {
                const double k  = kmin_th + ik * kstep_th;
                const double wk = simpsonWeight(ik, p_.nk);
                F2pi_k += wk * integrand(x, Q2, th, k);
            }
            F2pi_k *= (kstep_th / 3.0);

            F2piK += wx * F2pi_k;
        }
        F2piK *= (xstep_th / 3.0);

        // Accumulate first moment (Simpson composite over theta)
        F2pi0 += simpsonWeight(ith, p_.nth) * F2piK;

        ThetaScanPoint pt;
        pt.theta_deg   = th * 180.0 / PI;
        pt.F2piK       = F2piK;
        pt.F2piMoment  = 0.0;   // filled in second pass below
        out.thetaScan.push_back(pt);

        if (p_.verbose)
            std::cout << "  [theta-scan]  theta=" << pt.theta_deg
                      << " deg   F2piK=" << F2piK << "\n";
    }

    out.F2pi0_theta = (thstep / 3.0) * F2pi0;

    // Running moment using composite Simpson (consistent with first moment)
    {
        int np = static_cast<int>(out.thetaScan.size());
        double cumulative = 0.0;
        for (int i = 0; i < np; ++i) {
            cumulative += simpsonWeight(i, np - 1) * out.thetaScan[i].F2piK;
            out.thetaScan[i].F2piMoment = (thstep / 3.0) * cumulative;
        }
    }

    if (p_.verbose)
        std::cout << "  [theta-scan]  First moment of F2piK(theta_h) = "
                  << out.F2pi0_theta << "\n";
}

// ----------------------------------------------------------
//  runAll -- run scans selected by p_.scanMode
// ----------------------------------------------------------
ScanResults Calculator::runAll() const
{
    ScanResults res;
    if (p_.scanMode == "x"     || p_.scanMode == "all") {
        if (p_.verbose) std::cout << "\n=== Running x-scan ===\n";
        runXScan(res);
    }
    if (p_.scanMode == "theta" || p_.scanMode == "all") {
        if (p_.verbose) std::cout << "\n=== Running theta-scan ===\n";
        runThetaScan(res);
    }
    if (p_.scanMode == "kbin"  || p_.scanMode == "all") {
        if (p_.verbose) std::cout << "\n=== Running kbin-scan ===\n";
        runKBinScan(res);
    }
    return res;
}

// ----------------------------------------------------------
//  eval -- single-point F2pi evaluation for generator use
//
//  Runs the kbin (y,kT) double integral at a single x point.
//  If Q2 <= 0, it is computed from (E, x, theta_e) internally.
//  All RunParams kinematic cuts and GRV/CTEQ6 settings apply.
//  Output is always silent (verbose flag is not used here).
// ----------------------------------------------------------
double Calculator::eval(double x, double Q2) const
{
    if (Q2 <= 0.0)
        Q2 = calcQ2(p_.E, x, p_.theta_e);
    if (Q2 <= 0.0 || x <= 0.0 || x >= 1.0)
        return 0.0;

    const double kTmax  = 10.0;
    const double kTstep = kTmax / p_.nkT;
    const double y0     = x;
    const double ystep  = (1.0 - y0) / p_.ny;

    double F2piK = 0.0;

    for (int iy = 0; iy <= p_.ny; ++iy) {
        const double y  = y0 + iy * ystep;
        const double wy = simpsonWeight(iy, p_.ny);

        double fpi_intk = 0.0;
        for (int ikT = 0; ikT <= p_.nkT; ++ikT) {
            const double kT  = ikT * kTstep;
            const double wkT = simpsonWeight(ikT, p_.nkT);

            const double km  = kMag  (kT, y);
            const double cph = cosPhi(kT, y);
            if (cph > p_.cosph_max) continue;
            if (cph < p_.cosph_min) continue;
            if (km  < p_.kmin)      continue;
            if (km  > p_.kmax)      continue;

            fpi_intk += wkT * 2.0 * kT * splitting(y, kT);
        }
        fpi_intk *= (kTstep / 3.0);

        if (y > 0.0 && x / y < 1.0)
            F2piK += wy * fpi_intk * f2pi(x / y, Q2);
    }
    return F2piK * (ystep / 3.0);
}


//
//  Direct translation of the Fortran expression (line 153-154):
//    kmag = sqrt( kT² + [kT² + (1-(1-y)²)mN²]² / [4mN²(1-y)²] )
//
//  Physical meaning: given the transverse momentum kT and the
//  meson IMF fraction y, reconstruct the full |k| of the recoil
//  nucleon using 4-momentum conservation.
// ----------------------------------------------------------
double Calculator::kMag(double kT, double y)
{
    const double omy  = 1.0 - y;
    const double num  = kT*kT + (1.0 - omy*omy) * mN*mN;
    return std::sqrt( kT*kT + num*num / (4.0 * mN*mN * omy*omy) );
}

// ----------------------------------------------------------
//  cosPhi -- proton polar angle cosine from (kT, y)
//
//  Direct translation of the Fortran expression (lines 156-158):
//    cosph = [kT² + (1-(1-y)²)mN²] / [2mN(1-y) * kmag]
// ----------------------------------------------------------
double Calculator::cosPhi(double kT, double y)
{
    const double omy = 1.0 - y;
    const double km  = kMag(kT, y);
    if (km <= 0.0) return 0.0;
    return (kT*kT + (1.0 - omy*omy)*mN*mN) / (2.0 * mN * omy * km);
}

// ----------------------------------------------------------
//  runKBinScan -- mirrors f2pi_sub.f90
//
//  Output loop : x  [xmin, xmax]             (nx+1 points)
//  Outer loop  : y  [x, 1]                   (ny steps, Simpson)
//  Inner loop  : kT [0, kTmax]               (nkT steps, Simpson)
//
//  At each (y, kT):
//   1. Compute kmag(kT,y) and cosph(kT,y)
//   2. Reject if cosph outside [cosph_min, cosph_max]  (angle cut)
//   3. Reject if kmag outside [kmin, kmax]              (|k| bin)
//   4. Accumulate 2*kT * splitFunc(y,kT) into kT integral
//  Then:
//   5. Multiply kT integral by F2piGRV(x/y, Q2) → y integrand
//   6. Accumulate y integral → F2piK(x)
//
//  This exactly reproduces the Fortran f2pi_sub calculation.
//  kT is integrated from 0 to kTmax=10 GeV (Fortran default).
// ----------------------------------------------------------
void Calculator::runKBinScan(ScanResults &out) const
{
    out.xScan.clear();
    out.xScan.reserve(p_.nx + 1);

    const double xstep  = (p_.xmax - p_.xmin) / p_.nx;
    const double kTmax  = 10.0;                           // Fortran: kTmaxk = 10
    const double kTstep = kTmax / p_.nkT;

    double F2pi0 = 0.0;

    if (p_.verbose)
        std::cout << std::scientific << std::setprecision(5);

    for (int ix = 0; ix <= p_.nx; ++ix) {
        const double x  = p_.xmin + ix * xstep;
        const double Q2 = calcQ2(p_.E, x, p_.theta_e);

        // Small-x warning: at theta_e=12 deg, Q2 can be negative for x<0.01
        if (Q2 <= 0.0 && p_.verbose)
            std::cerr << "  [kbin-scan] WARNING: Q2<=0 at x=" << x
                      << " (theta_e=" << p_.theta_e << " deg) -- point skipped.\n";

        const double y0    = x;
        const double y1    = 1.0;
        const double ystep = (y1 - y0) / p_.ny;

        double F2piK = 0.0;

        for (int iy = 0; iy <= p_.ny; ++iy) {
            const double y  = y0 + iy * ystep;
            const double wy = simpsonWeight(iy, p_.ny);

            double fpi_intk = 0.0;

            for (int ikT = 0; ikT <= p_.nkT; ++ikT) {
                const double kT  = ikT * kTstep;
                const double wkT = simpsonWeight(ikT, p_.nkT);

                const double km  = kMag  (kT, y);
                const double cph = cosPhi(kT, y);

                if (cph > p_.cosph_max) continue;
                if (cph < p_.cosph_min) continue;
                if (km  < p_.kmin)      continue;
                if (km  > p_.kmax)      continue;

                fpi_intk += wkT * 2.0 * kT * splitting(y, kT);
            }
            fpi_intk *= (kTstep / 3.0);

            // pion structure function at x_pi = x/y (set chosen in f2pi())
            if (y > 0.0 && x / y < 1.0 && Q2 > 0.0)
                F2piK += wy * fpi_intk * f2pi(x / y, Q2);
        }
        F2piK *= (ystep / 3.0);

        F2pi0 += simpsonWeight(ix, p_.nx) * F2piK;

        XScanPoint pt;
        pt.x          = x;
        pt.Q2         = Q2;
        pt.F2piK      = F2piK;
        pt.F2piMoment = 0.0;
        pt.F2n        = 0.0;
        pt.ratio      = 0.0;

#ifdef HAVE_CTEQ6
        if (Q2 > 0.0) {
            pt.F2n = (p_.nucleon == NUCLEON_PROTON)
                   ? CTEQ6::F2proton (x, Q2, false)
                   : CTEQ6::F2neutron(x, Q2, false);
            if (pt.F2n > 0.0) pt.ratio = pt.F2piK / pt.F2n;
        }
#endif
        out.xScan.push_back(pt);

        if (p_.verbose)
            std::cout << "  [kbin-scan]  x=" << x
                      << "  Q2=" << Q2
                      << "  F2piK=" << F2piK
#ifdef HAVE_CTEQ6
                      << "  ratio=" << pt.ratio
#endif
                      << "\n";
    }

    out.F2pi0_x = (xstep / 3.0) * F2pi0;

    {
        int np = static_cast<int>(out.xScan.size());
        double cumulative = 0.0;
        for (int i = 0; i < np; ++i) {
            cumulative += simpsonWeight(i, np - 1) * out.xScan[i].F2piK;
            out.xScan[i].F2piMoment = (xstep / 3.0) * cumulative;
        }
    }

    if (p_.verbose)
        std::cout << "  [kbin-scan]  First moment = " << out.F2pi0_x
                  << "  |k| bin=[" << p_.kmin << "," << p_.kmax << "] GeV\n";
}

} // namespace PionCloud
