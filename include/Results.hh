#pragma once
// ============================================================
//  Results.hh
//  Pion Cloud Model -- data containers for scan results
//
//  These plain structs are the only interface between the
//  Calculator (physics) and the OutputWriter (I/O).  They
//  carry no ROOT or file-system dependencies and can be used
//  directly in embedded / library contexts.
// ============================================================

#include <vector>

namespace PionCloud {

// ----------------------------------------------------------
//  XScanPoint
//  One output point from runXScan (outer loop = x).
// ----------------------------------------------------------
struct XScanPoint {
    double x;             // Bjorken-x
    double Q2;            // Q² at this x [GeV²]
    double F2piK;         // F₂^π(x) integrated over (θ_h, |k|)
    double F2piMoment;    // running ∫_{xmin}^{x} F₂^π(x') dx'  (trapezoidal)
    double F2n;           // F₂ⁿ(x,Q²) from CTEQ6 + charge symmetry
                          //   = 0 when CTEQ6 is not available
    double ratio;         // F₂^π / F₂ⁿ  (0 when F2n=0 or CTEQ6 absent)
};

// ----------------------------------------------------------
//  ThetaScanPoint
//  One output point from runThetaScan (outer loop = θ_h).
// ----------------------------------------------------------
struct ThetaScanPoint {
    double theta_deg;     // hadron production angle [degrees]
    double F2piK;         // F₂^π(θ_h) integrated over (x, |k|)
    double F2piMoment;    // running ∫_{θmin}^{θ} F₂^π(θ') dθ'  (trapezoidal)
};

// ----------------------------------------------------------
//  ScanResults
//  Output container populated by Calculator::runXScan /
//  runThetaScan / runAll.  Fields are empty / zero until
//  the corresponding scan has been run.
// ----------------------------------------------------------
struct ScanResults {
    std::vector<XScanPoint>     xScan;
    std::vector<ThetaScanPoint> thetaScan;

    // First moments computed with composite Simpson's rule,
    // matching the Fortran PRINT* statements exactly.
    double F2pi0_x     = 0.0;   // ∫ F₂^π(x) dx
    double F2pi0_theta = 0.0;   // ∫ F₂^π(θ_h) dθ_h
};

} // namespace PionCloud
