#pragma once
// ============================================================
//  Calculator.hh
//  Pion Cloud Model -- main calculation engine
//
//  Three scan modes, each a direct C++ port of a Fortran program:
//
//    runXScan()     -- mirrors 3Var_x.f
//                      outer: x, integrated: (theta_h, |k|)
//                      → F₂^π(x) for a given |k| tagging window
//
//    runThetaScan() -- mirrors 3Var_theta.f
//                      outer: theta_h, integrated: (x, |k|)
//                      → F₂^π(theta_h) for Fig.47 style plots
//
//    runKBinScan()  -- mirrors f2pi_sub.f90
//                      outer: x, integrated: (y, kT)
//                      applies kmag window [kmin,kmax] and
//                      proton angle acceptance [cosph_min,cosph_max]
//                      → F₂^π(x) for the four k-bin curves
//
//  All three store results into ScanResults and share the same
//  RunParams.  scan_mode = "x" | "theta" | "kbin" | "all"
//  selects which ones run.  "all" runs all three.
//
//  The class has no I/O; output is handled by OutputWriter.
// ============================================================

#include "PhysicsParams.hh"
#include "Results.hh"

namespace PionCloud {

class Calculator {
public:
    explicit Calculator(const RunParams &p);

    // ---- individual scans ----
    void runXScan    (ScanResults &out) const;
    void runThetaScan(ScanResults &out) const;
    void runKBinScan (ScanResults &out) const;

    // ---- convenience: run scans selected by p_.scanMode ----
    ScanResults runAll() const;

    // ---- single-point evaluation (generator/library interface) ----
    //
    // eval(x, Q2) returns F₂^π at one (x, Q²) point using the kbin
    // integration strategy (y, kT) with the cuts in RunParams.
    // Equivalent to running runKBinScan at a single x point with
    // verbose=false, without modifying any stored state.
    //
    // This is the primary entry point for embedded use:
    //   Calculator calc(p);
    //   double f2pi = calc.eval(0.1, 2.0);
    //
    // Q2 is accepted as a parameter so callers that already know Q2
    // (e.g. from their own kinematics) can bypass the internal
    // calcQ2() call.  Pass Q2 <= 0 to let the Calculator compute it
    // from (E, x, theta_e) via calcQ2().
    double eval(double x, double Q2 = -1.0) const;

private:
    RunParams p_;

    // (theta_h, k) integrand for runXScan / runThetaScan
    double integrand(double x, double Q2,
                     double theta_h, double k) const;

    // Compute 3D momentum magnitude from (kT, y) -- f2pi_sub formula
    static double kMag  (double kT, double y);
    // Compute proton polar angle cosine from (kT, y) -- f2pi_sub formula
    static double cosPhi(double kT, double y);
};

} // namespace PionCloud
