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
//
//  ------------------------------------------------------------
//  HOW TO MODIFY  (v1.4)
//  ------------------------------------------------------------
//  Every scan evaluates the same physics,
//      dF2 = f_MB(y, kT) * F2^M(x/y, Q2)   [* Jacobian for (|k|,theta)],
//  integrated over a different region / set of variables.
//
//  * New splitting function or meson structure function:
//      do NOT edit the scans.  Implement SplittingModel / PionSFModel
//      (Models.hh) and inject it (2nd constructor) or register it by
//      name (input keys splitting_model / pion_sf_model).
//      All scans, output and the C interface then use it.
//  * New built-in channel: SplittingFunctions.hh ("HOW TO ADD A CHANNEL").
//  * Different Q2 prescription: calcQ2() (Kinematics.hh) is called in
//      runXScan, runThetaScan, runKBinScan and eval; or pass Q2 to eval().
//  * New cut / acceptance in the kbin scan: the `continue` tests on
//      kmag and cosph in runKBinScan() and eval() (keep them identical).
//  * Integration grids: all loops use composite Simpson via
//      simpsonWeight() (Kinematics.hh), valid for any number of steps.
//  * Which RunParams fields each scan reads:
//      runXScan     : xmin,xmax,nx | thmin_deg,thmax_deg,nth_int | kmin,kmax,nk
//      runThetaScan : thmin_deg,thmax_deg,nth (OUTPUT grid) |
//                     xmin_th,xmax_th,nx_int | kmin_th,kmax_th,nk
//      runKBinScan  : xmin,xmax,nx | ny | nkT (kT in [0,10] GeV) |
//                     kmin,kmax (|k| window) | cosph_min,cosph_max
//      all          : E, theta_e (-> Q2), model fields (flag, L/Ld, typ, dis,
//                     grv_version, pion_strange, *_model)
//    NOTE: theta_min/max are an OUTPUT range in the theta scan but an
//    INTEGRATION range in the x scan; kmin/kmax are the x-scan |k| range
//    and the kbin window; nk is shared by the x and theta scans.
//  * After any change: ctest, then docs/validation/run_compare.sh and
//    compare with the published curves (see CONTRIBUTING.md).
// ============================================================

#include "PhysicsParams.hh"
#include "Results.hh"
#include "Models.hh"
#include "Quadrature.hh"

#include <memory>

namespace PionCloud {

class Calculator {
public:
    explicit Calculator(const RunParams &p);

    // Inject user models (v1.4).  A null pointer means "use the model
    // selected by p" (built-in or registered by name, see Models.hh).
    Calculator(const RunParams &p,
               std::shared_ptr<const SplittingModel> split,
               std::shared_ptr<const PionSFModel>    pionSF = nullptr);

    const SplittingModel &splittingModel() const { return *split_; }
    const PionSFModel    &pionSFModel()    const { return *sf_; }

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
    std::shared_ptr<const SplittingModel> split_;   // f_MB(y, kT)
    std::shared_ptr<const PionSFModel>    sf_;      // F2^M(x_M, Q2)
    GaussLegendre                         gl_;      // nodes for KBIN_EXACT
    double chargeFactor_ = 1.0;   // isospin sum factor (total_charge_states)

    // Q2 at Bjorken x: q2_fixed if > 0, else calcQ2(E, x, theta_e)
    double Q2at(double x) const;

    // kbin quantity F2^(MB)(x; |k| window, angle window) at one (x, Q2):
    double kbinF2     (double x, double Q2) const;  // dispatches on kbin_method
    double kbinF2Grid (double x, double Q2) const;  // Fortran f2_pi_sub algorithm
    double kbinF2Exact(double x, double Q2) const;  // exact region + Gauss-Legendre

    // ---- physics ingredients (one place each; v1.4) ----
    // Meson-baryon splitting function f_MB(y, kT) of the selected channel
    // (p_.flag), with the cut-off that channel uses (L or Ld).
    double splitting(double y, double kT) const;
    // Pion structure function F2pi(x_pi, Q2) of the selected PDF set
    // (p_.grv_version), with/without strange quarks (p_.pion_strange).
    double f2pi(double xpi, double Q2) const;

    // (theta_h, k) integrand for runXScan / runThetaScan
    double integrand(double x, double Q2,
                     double theta_h, double k) const;

    // Compute 3D momentum magnitude from (kT, y) -- f2pi_sub formula
    static double kMag  (double kT, double y);
    // Compute proton polar angle cosine from (kT, y) -- f2pi_sub formula
    static double cosPhi(double kT, double y);
};

} // namespace PionCloud
