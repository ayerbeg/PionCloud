#pragma once
// ============================================================
//  PhysicsParams.hh
//  Pion Cloud Model -- Physics constants and run parameters
//
//  All masses in GeV, angles in radians (converted at input).
// ============================================================

#include <cmath>
#include <string>

namespace PionCloud {

// ----------------------------------------------------------
//  Mathematical constant
// ----------------------------------------------------------
constexpr double PI = M_PI;

// ----------------------------------------------------------
//  Physical masses  [GeV]
// ----------------------------------------------------------
constexpr double mN    = 0.93891897;   // nucleon mass
constexpr double mpi   = 0.13957018;   // pion mass
constexpr double mrho  = 0.7754;       // rho meson
constexpr double mDelta= 1.232;        // Delta baryon

// ----------------------------------------------------------
//  Form-factor type codes
// ----------------------------------------------------------
enum FormFactorType {
    FF_MONOPOLE    = 0,
    FF_DIPOLE      = 1,
    FF_EXP_S       = 2,   // s-dependent exponential
    FF_COV_DIPOLE  = 3,
    FF_DIPOLE_S_CH = 4,   // dipole, s-channel Lambda exchange
    FF_EXP_T       = 5,   // t-dependent exponential
    FF_PAULI_VILLAR= 6
};

// ----------------------------------------------------------
//  Splitting-function / meson-baryon FLAG codes
// ----------------------------------------------------------
enum ContributionFlag {
    FLAG_PION     = 0,   // pi-N  (J = 0+1/2)
    FLAG_RHO      = 1,   // rho-N (J = 1+1/2)
    FLAG_PI_DELTA = 2    // pi-Delta (J = 0+3/2)
};

// ----------------------------------------------------------
//  Dissociation channel codes
// ----------------------------------------------------------
enum DisChannel {
    DIS_CHARGE_EXCHANGE = 0,   // N -> pi^- + p  (charge exchange)
    DIS_NEUTRAL         = 1    // N -> pi^0 + N  (neutral)
};

// ----------------------------------------------------------
//  Nucleon type for F2 denominator in ratio output
// ----------------------------------------------------------
enum NucleonType {
    NUCLEON_NEUTRON = 0,   // F2n via charge symmetry (default, matches Fortran)
    NUCLEON_PROTON  = 1    // F2p directly from CTEQ6
};

// ----------------------------------------------------------
//  GRV pion PDF version
// ----------------------------------------------------------
enum GRVVersion {
    GRV_92 = 92,   // GRV92 LO, Z.Phys.C53(1992)651 (default)
    GRV_99 = 99    // GRV99 LO update (in f2_pi_sub.f90 as subroutine GRV99)
};

// ----------------------------------------------------------
//  Precision preset for kbin scan
//  Controls ny and nkT together.
// ----------------------------------------------------------
enum PrecisionPreset {
    PRECISION_CUSTOM      = 0,   // use ny/nkT directly (default)
    PRECISION_FAST        = 1,   // ny=20,  nkT=500   -- quick checks
    PRECISION_NORMAL      = 2,   // ny=50,  nkT=2000  -- development
    PRECISION_PUBLICATION = 3    // ny=100, nkT=50000 -- matches Fortran
};
enum OutputMode {
    OUT_NONE  = 0,
    OUT_ASCII = 1 << 0,   // write ASCII summary + data tables
    OUT_ROOT  = 1 << 1,   // write ROOT file with TGraphs/TH1
    OUT_ALL   = OUT_ASCII | OUT_ROOT
};

// ----------------------------------------------------------
//  Run-parameters struct
//  Mirrors the hard-coded values in the original Fortran
//  programs but collected here for easy modification.
// ----------------------------------------------------------
struct RunParams {
    // ---- beam/kinematics ----
    double E        = 11.0;    // electron beam energy [GeV]
    double theta_e  = 35.0;    // electron scattering angle [deg]

    // ---- lambda cut-offs ----
    double L        = 1.33;    // renorm. cut-off for pi/rho-N vertex [GeV]
    double Ld       = 1.39;    // renorm. cut-off for Delta vertex   [GeV]

    // ---- form factor & channel choices ----
    FormFactorType  typ  = FF_PAULI_VILLAR;
    ContributionFlag flag = FLAG_PION;
    DisChannel       dis  = DIS_CHARGE_EXCHANGE;

    // ---- hadron kinematics (from subroutine / f2_pi_sub usage) ----
    double pH         = 0.325;  // produced hadron momentum [GeV]
    double alpha1_deg = 30.0;   // hadron prod. angle lower bound [deg]
    double alpha2_deg = 70.0;   // hadron prod. angle upper bound [deg]

    // ---- x-scan (3Var_x) ----
    double xmin     = 0.02;
    double xmax     = 0.30;
    int    nx       = 100;       // number of output x points
    int    nth_int  = 100;       // internal theta_h integration steps (x-scan)

    // ---- theta_h-scan (3Var_theta) ----
    double thmin_deg  = 0.0;    // hadron production angle lower bound [deg]
    double thmax_deg  = 100.0;  // hadron production angle upper bound [deg]
    int    nth        = 100;    // number of output theta points
    int    nx_int     = 100;    // internal x integration steps (theta-scan)
    // x integration bounds for theta-scan
    // Default: [~0, 0.6] -- full accessible range, matching 3Var_theta.f
    // Right panel of Fig.47: set xmin_th=0.05, xmax_th=0.6
    double xmin_th    = 1.0e-4; // avoids x=0 singularity in Q2
    double xmax_th    = 0.6;
    // |k| integration bounds for theta-scan
    // Default: [~0, 0.5 GeV] -- full range, matching 3Var_theta.f
    // Right panel of Fig.47: set kmin_th=0.060, kmax_th=0.250
    double kmin_th    = 1.0e-4; // avoids k=0 (rho_PS vanishes there)
    double kmax_th    = 0.5;

    // ---- |k| integration (shared) ----
    double kmin     = 0.050;   // [GeV] -- x-scan theta window / kbin window
    double kmax     = 0.100;   // [GeV]
    int    nk       = 100;

    // ---- kbin-scan (f2pi_sub mode) ----
    // Integrates over (y, kT) at fixed x, applying a |k| window and
    // a proton polar angle acceptance cut on cosph = cos(phi_k).
    // Reproduces the four k-bin curves from f2pi_sub / f2_pi_ranges.c.
    //
    // cosph_min = cos(70°) = 0.342  (70° upper angle cut)
    // cosph_max = cos(30°) = 0.866  (30° lower angle cut)
    // These are the values hardcoded in f2pi_sub.f90 lines 170,174.
    // Exposing them here lets you change the angular acceptance window.
    double cosph_min  = 0.342;   // cos(70°): reject cosph < this
    double cosph_max  = 0.866;   // cos(30°): reject cosph > this
    int    ny         = 100;     // steps in y integration (f2pi_sub used 100)
    int    nkT        = 50000;   // steps in kT integration (f2pi_sub used 50000)
                                 // reduce to ~1000 for fast tests, 50000 for production

    // ---- nucleon type for F2 denominator ----
    // NUCLEON_NEUTRON (default): F2n via charge symmetry -- matches Fortran
    // NUCLEON_PROTON:            F2p directly from CTEQ6
    NucleonType nucleon  = NUCLEON_NEUTRON;

    // ---- GRV pion PDF version ----
    // GRV_92 (default): GRV92 LO, matches original Fortran 3Var_x.f
    // GRV_99:           GRV99 LO update from f2_pi_sub.f90
    GRVVersion  grv_version = GRV_92;

    // ---- kbin precision preset ----
    // PRECISION_CUSTOM (default): use ny/nkT fields directly
    // PRECISION_FAST:             ny=20,  nkT=500
    // PRECISION_NORMAL:           ny=50,  nkT=2000
    // PRECISION_PUBLICATION:      ny=100, nkT=50000 (matches Fortran)
    PrecisionPreset precision = PRECISION_CUSTOM;

    // ---- inclusive mode shortcut ----
    // When true, overrides cosph_min=0, cosph_max=1, kmin=0.001, kmax=10.0
    // Reproduces the fully inclusive F2^(piN)(x) curve (no acceptance cuts).
    bool inclusive_mode = false;

    // ---- verbosity ----
    // true (default): print per-point progress to stdout
    // false:          suppress all stdout (for library/generator embedding)
    bool verbose = true;

    // ---- CTEQ6 table directory (generator/library interface) ----
    // SetCtq6 opens the .tbl file by name in the current working directory.
    // If cteq6_table_dir is non-empty, PionCloud chdir()s to that
    // directory before calling CTEQ6::init() and back afterwards.
    // Example: p.cteq6_table_dir = "/data/cteq6_tables"
    // Leave empty (default) to use the current working directory.
    std::string cteq6_table_dir = "";

    // ---- run control ----
    std::string scanMode    = "all"; // "x", "theta", "kbin", or "all"

    // ---- output ----
    int         outputMode  = OUT_ALL;
    std::string outBaseName = "PionCloud";
};

} // namespace PionCloud
