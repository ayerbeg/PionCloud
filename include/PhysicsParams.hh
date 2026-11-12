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
    FLAG_PION      = 0,   // pi-N      (J = 0+1/2)
    FLAG_RHO       = 1,   // rho-N     (J = 1+1/2)
    FLAG_PI_DELTA  = 2,   // pi-Delta  (J = 0+3/2)   uses Ld
    FLAG_RHO_DELTA = 3    // rho-Delta (J = 1+3/2)   uses Ld   (selectable since v1.4;
                          //   f_RhoDel existed in v1.3 but could not be chosen)
};

// ----------------------------------------------------------
//  Dissociation channel codes
// ----------------------------------------------------------
// Meaning per channel (isospin factor applied in SplittingFunctions.hh):
//   pi-N, rho-N : 0 -> n -> M^- p  (or p -> M^+ n), factor 2
//                 1 -> N -> M^0 N,                   factor 1
//   pi-Delta, rho-Delta : only 0 is implemented, meaning the SUM of the
//                 charge states of the Fortran (factor 2/3 x coupling);
//                 1 has zero coupling in the Fortran and is REJECTED since
//                 v1.4 (it silently returned F2 = 0 before).
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
    GRV_92 = 92,   // GRV92 LO pion, Z.Phys.C53(1992)651 (default, = Fortran)
    GRV_99 = 99    // GRS99 LO pion, Gluck-Reya-Schienbein EPJ C10(1999)313
                   // (corrected in v1.4 -- see GRV.hh and CHANGELOG.md)
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
// ----------------------------------------------------------
//  kbin integration method (v1.4)
// ----------------------------------------------------------
enum KBinMethod {
    KBIN_GRID  = 0,   // ORIGINAL f2_pi_sub.f90 algorithm: uniform kT grid on
                      // [0,10] GeV, ny x nkT Simpson, cuts applied by rejecting
                      // grid points.  Select with kbin_method = grid / --kbin-grid.
    KBIN_EXACT = 1    // NEW in v1.4, DEFAULT: exact |k|/angle region in
                      // u = kT^2 + Gauss-Legendre
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
    // Defaults (v1.4): s-dependent exponential form factor with the HSS
    // central values  Lambda_piN = Lambda_rhoN = 1.56 GeV,
    // Lambda_piDelta = Lambda_rhoDelta = 1.39 GeV  (Holtmann, Szczurek, Speth,
    // NPA 596 (1996) 631, in the exp[(M^2-s)/Lambda^2] convention of the TDIS
    // proposal PR12-15-006 Eq. 14).  These reproduce the published TDIS figures.
    //
    // WHY CHANGED: v1.3 defaulted to Pauli-Villars with L = 1.33 GeV.  1.33 GeV
    // is the best fit of the s-dependent EXPONENTIAL in McKenney et al.,
    // PRD 93 (2016) 054011, Table I; the Pauli-Villars fit there is 0.27 GeV.
    // The pairing (inherited from 3Var_x.f) matches no published fit, Ld = 1.39
    // was already the HSS exp_s value, and the C interface default
    // (pioncloud_init_default) already used exp_s / 1.56.  Set form_factor and
    // lambda explicitly to get any other model.
    double L        = 1.56;    // cut-off for pi/rho-N vertex      [GeV]
    double Ld       = 1.39;    // cut-off for pi/rho-Delta vertex  [GeV]

    // ---- form factor & channel choices ----
    FormFactorType  typ  = FF_EXP_S;
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
    // GRV_99:           GRS99 LO (hep-ph/9903288); v1.3 had this set wrong
    //                   (valence x2, strange sea used as light sea) -- fixed v1.4
    GRVVersion  grv_version = GRV_92;

    // ---- strange quarks in the pion structure function ----
    // false (default, = Fortran): F2pi = 5/9 (xV + 2 xS)           (u, d only)
    // true:                       F2pi += 2/9 x sbar                (adds s, sbar)
    // GRV92 has an SU(3)-symmetric sea (sbar = S); GRS99 has its own sbar.
    bool pion_strange = false;

    // ---- Q2 prescription (v1.4) ----
    // q2_fixed <= 0 (default): Q2(x) from beam energy E and fixed electron
    //   angle theta_e (calcQ2, as in the Fortran).
    // q2_fixed  > 0: use this Q2 [GeV^2] at every x (model studies,
    //   comparison with theory at fixed scale).
    double q2_fixed = 0.0;

    // ---- charge states (v1.4) ----
    // false (default): one charge state, as selected by `dis` (Fortran).
    // true: sum over charge states for pi-N / rho-N
    //   (neutral x3, charge_exchange x3/2: isospin factors 1+2 = 3);
    //   pi-Delta / rho-Delta are already summed in the Fortran -> x1.
    //   This is the "fully integrated" F2^(piN) of the published figures.
    bool total_charge_states = false;

    // ---- user models (v1.4, see Models.hh) ----
    // Empty (default): built-in model selected by flag / grv_version.
    // Otherwise the name of a model registered with ModelRegistry
    // (input keys splitting_model / pion_sf_model).
    std::string splitting_model = "";
    std::string pion_sf_model   = "";

    // ---- kbin precision preset ----
    // PRECISION_CUSTOM (default): use ny/nkT fields directly
    // PRECISION_FAST:             ny=20,  nkT=500
    // PRECISION_NORMAL:           ny=50,  nkT=2000
    // PRECISION_PUBLICATION:      ny=100, nkT=50000 (matches Fortran)
    PrecisionPreset precision = PRECISION_CUSTOM;

    // ---- kbin integration method (NEW in v1.4) ----
    // KBIN_EXACT (DEFAULT since v1.4): solves the |k| window and angular
    //   cuts analytically in u = kT^2 at each y (at most two intervals) and
    //   integrates with Gauss-Legendre: no rejected points, no step-function
    //   integrand.  0.002% accuracy in ~0.3 s for a 50-point scan.
    //   gl_nodes nodes per panel; y in [x, ymax] split into gl_ypanels.
    // KBIN_GRID: the ORIGINAL Fortran algorithm (f2_pi_sub.f90), kept to
    //   reproduce Fortran / v1.3 results exactly.  Select it with
    //   kbin_method = grid  (input file),  --kbin-grid  (command line),
    //   pioncloud_set_kbin_method(0)  (C interface).
    //   ny / nkT / precision presets apply ONLY to this method.
    // Typical difference exact vs grid (Fortran settings ny=100, nkT=5000):
    //   < 1% in the bulk, several % close to the kinematic end point where
    //   the grid result is least accurate.
    KBinMethod kbin_method = KBIN_EXACT;
    int        gl_nodes    = 16;
    int        gl_ypanels  = 32;

    // ---- inclusive mode shortcut ----
    // When true, overrides cosph_min=0, cosph_max=1, kmin=0.001, kmax=10.0
    // Reproduces the fully inclusive F2^(piN)(x) curve of ONE charge state
    // (x3 for neutral, x3/2 for charge exchange gives the total).
    // "Inclusive" up to the |k| limits: recoils with |k| < 1 MeV or
    // |k| > 10 GeV are still excluded -- negligible in practice.
    // cosph >= 0 always holds (the recoil goes forward), so [0,1] is all angles.
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
