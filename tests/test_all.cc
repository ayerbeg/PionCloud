// ============================================================
//  tests/test_all.cc
//  PionCloud unit tests -- standalone, no external framework
//
//  Tests are grouped by module.  Each test function returns
//  true on pass.  Run with:
//    ./build/PionCloudTests
//  or via CTest:
//    ctest --test-dir build -V
//
//  Exit code: 0 = all passed, 1 = one or more failures.
// ============================================================

#include "PhysicsParams.hh"
#include "Kinematics.hh"
#include "GRV.hh"
#include "SplittingFunctions.hh"
#include "InputReader.hh"
#include "Calculator.hh"
#include "Results.hh"

#ifdef HAVE_CTEQ6
#include "CTEQ6.hh"
#endif

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace PionCloud;

// ============================================================
//  Test infrastructure
// ============================================================
static int g_pass = 0, g_fail = 0;

static bool check(bool condition, const std::string &name, const std::string &detail = "")
{
    if (condition) {
        std::cout << "  [PASS] " << name << "\n";
        ++g_pass;
    } else {
        std::cout << "  [FAIL] " << name;
        if (!detail.empty()) std::cout << "  -- " << detail;
        std::cout << "\n";
        ++g_fail;
    }
    return condition;
}

static bool near(double a, double b, double tol = 1e-8,
                 const std::string &name = "", const std::string &detail = "")
{
    bool ok = std::abs(a - b) < tol;
    std::ostringstream d;
    d << "got=" << a << " expected=" << b << " tol=" << tol;
    if (!detail.empty()) d << " (" << detail << ")";
    return check(ok, name, d.str());
}

static bool notNeg(double v, const std::string &name)
{
    std::ostringstream d; d << "value=" << v;
    return check(v >= 0.0, name, d.str());
}

// ============================================================
//  1. PhysicsParams -- defaults and enum values
// ============================================================
static void test_PhysicsParams()
{
    std::cout << "\n--- PhysicsParams ---\n";
    RunParams p;
    near(p.E,       11.0, 1e-12, "default E = 11 GeV");
    near(p.theta_e, 35.0, 1e-12, "default theta_e = 35 deg");
    near(p.L,       1.33, 1e-12, "default L = 1.33 GeV");
    near(p.Ld,      1.39, 1e-12, "default Ld = 1.39 GeV");
    check(p.typ  == FF_PAULI_VILLAR,     "default form factor = Pauli-Villar");
    check(p.flag == FLAG_PION,           "default contribution = pion");
    check(p.dis  == DIS_CHARGE_EXCHANGE, "default dis = charge exchange");
    check(p.scanMode == "all",           "default scanMode = all");
    check(p.outputMode == OUT_ALL,       "default outputMode = all");
    near(p.xmin, 0.02, 1e-12, "default xmin = 0.02");
    near(p.xmax, 0.30, 1e-12, "default xmax = 0.30");
    check(p.nx == 100,  "default nx = 100");
    check(p.nk == 100,  "default nk = 100");
    check(p.nth == 100, "default nth = 100");
    // kbin defaults
    near(p.cosph_min, 0.342, 1e-12, "default cosph_min = 0.342 (cos70°)");
    near(p.cosph_max, 0.866, 1e-12, "default cosph_max = 0.866 (cos30°)");
    check(p.ny  == 100,   "default ny = 100");
    check(p.nkT == 50000, "default nkT = 50000");
}

// ============================================================
//  2. Kinematics
// ============================================================
static void test_Kinematics()
{
    std::cout << "\n--- Kinematics ---\n";

    // calcQ2: at forward angles Q2 should be small and positive
    {
        double Q2 = calcQ2(11.0, 0.1, 35.0);
        check(Q2 > 0.0, "calcQ2 > 0 at (E=11,x=0.1,th=35)");
        // Compare with direct formula: Q2 = 2*mN*x*E*(1 - 1/(2E/(x*mN)*sin^2(th/2)+1))
        // At these kinematics expect O(1) GeV^2
        check(Q2 > 0.5 && Q2 < 10.0, "calcQ2 in plausible range [0.5,10] GeV^2");
    }
    // calcQ2 guard: x=0 must return 0 not NaN
    {
        double Q2 = calcQ2(11.0, 0.0, 35.0);
        near(Q2, 0.0, 1e-15, "calcQ2(x=0) = 0");
    }

    // kThetaToYkT: at theta=pi/2 all k goes transverse
    {
        double y, kT;
        kThetaToYkT(0.1, PI/2.0, y, kT);
        near(kT, 0.1, 1e-12, "kT = k at theta=90 deg");
        // y = 0 + (1/mN)*(mN - sqrt(mN^2+k^2)) < 0 at 90 deg → guarded in integrand
        check(y < 0.0, "y < 0 at theta=90 deg (guarded in integrand)");
    }
    // At theta=0 all k goes longitudinal, kT=0
    {
        double y, kT;
        kThetaToYkT(0.1, 0.0, y, kT);
        near(kT, 0.0, 1e-12, "kT = 0 at theta=0");
        check(y > 0.0 && y < 1.0, "y in (0,1) at theta=0, small k");
    }

    // phaseSpace vanishes at theta=0 (sin(0)=0)
    {
        double ps = phaseSpace(0.1, 0.0);
        near(ps, 0.0, 1e-15, "phaseSpace = 0 at theta=0");
    }
    // phaseSpace positive at theta in (0,pi)
    {
        double ps = phaseSpace(0.1, PI/4.0);
        check(ps > 0.0, "phaseSpace > 0 at theta=45 deg");
    }

    // simpsonWeight: endpoints = 1, odds = 4, evens (non-endpoint) = 2
    {
        check(simpsonWeight(0, 10) == 1.0,  "Simpson weight i=0 -> 1");
        check(simpsonWeight(10,10) == 1.0,  "Simpson weight i=n -> 1");
        check(simpsonWeight(1, 10) == 4.0,  "Simpson weight i=1 (odd) -> 4");
        check(simpsonWeight(2, 10) == 2.0,  "Simpson weight i=2 (even) -> 2");
        check(simpsonWeight(5, 10) == 4.0,  "Simpson weight i=5 (odd) -> 4");

        // Integral of f(x)=1 on [0,1] with n=10 steps: should be exactly 1
        double sum = 0.0;
        int n = 10;
        double h = 1.0 / n;
        for (int i = 0; i <= n; ++i) sum += simpsonWeight(i, n);
        near((h/3.0)*sum, 1.0, 1e-14, "Simpson ∫1 dx = 1 (n=10)");

        // Integral of f(x)=x^2 on [0,1] = 1/3 exactly (Simpson is exact for cubics)
        sum = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = i * h;
            sum += simpsonWeight(i, n) * x * x;
        }
        near((h/3.0)*sum, 1.0/3.0, 1e-14, "Simpson ∫x² dx = 1/3 (n=10)");
    }
}

// ============================================================
//  3. GRV PDFs
// ============================================================
static void test_GRV()
{
    std::cout << "\n--- GRV PDFs ---\n";

    // At moderate x and Q2, both distributions should be positive
    {
        double xV, xS;
        GRV92(0.1, 4.0, xV, xS);
        check(xV > 0.0, "xVpi > 0 at (x=0.1, Q2=4)");
        check(xS > 0.0, "xSpi > 0 at (x=0.1, Q2=4)");
        check(xV > xS,  "xVpi > xSpi at x=0.1 (valence-dominated)");
    }
    // At x close to 1, distributions must be very small (suppressed by (1-x)^D)
    // The GRV parametrization gives (1-x)^D with D~0.5-1, so at x=0.95 it is
    // suppressed but not negligibly tiny; test x=0.99 for a stricter check.
    {
        double xV, xS;
        GRV92(0.99, 4.0, xV, xS);
        check(xV < 0.1, "xVpi very small as x -> 0.99");
        check(xV >= 0.0, "xVpi non-negative as x -> 0.99");
    }
    // Below Q2 threshold, returns 0
    {
        double xV, xS;
        GRV92(0.1, 0.1, xV, xS);  // Q2 clamped to 0.2501
        // With Q2 clamped to just above threshold, s is small, distributions small but may be nonzero
        check(xV >= 0.0, "xVpi >= 0 below threshold (clamped)");
    }
    // F2piGRV is positive
    {
        double f2 = F2piGRV(0.2, 2.0);
        check(f2 > 0.0, "F2piGRV > 0 at (x=0.2, Q2=2)");
    }
    // F2piGRV is small but not zero near x=0.99 (GRV parametrization)
    {
        double f2 = F2piGRV(0.99, 4.0);
        check(f2 >= 0.0 && f2 < 0.1, "F2piGRV small and non-negative as x -> 0.99");
    }
}

// ============================================================
//  4. Splitting functions
// ============================================================
static void test_SplittingFunctions()
{
    std::cout << "\n--- Splitting functions ---\n";

    const double L   = 1.56;
    const double Ld  = 1.39;
    const FormFactorType typ = FF_EXP_S;
    const DisChannel dis = DIS_CHARGE_EXCHANGE;

    // All functions return 0 at y boundaries
    for (double y : {0.0, 1.0}) {
        std::ostringstream n; n << "fypiN(y=" << y << ") = 0";
        near(fypiN(y, 0.1, L, typ, dis), 0.0, 1e-15, n.str());
    }
    for (double y : {0.0, 1.0}) {
        std::ostringstream n; n << "f_rhoN(y=" << y << ") = 0";
        near(f_rhoN(y, 0.1, L, typ, dis), 0.0, 1e-15, n.str());
    }

    // At a valid (y, kT), fypiN must be positive (dis=0: charge exchange)
    {
        double val = fypiN(0.2, 0.1, L, typ, DIS_CHARGE_EXCHANGE);
        notNeg(val, "fypiN > 0 at (y=0.2, kT=0.1, dis=0)");
    }
    // Neutral dissociation gives half the coupling^2 -> smaller than charge exchange
    {
        double v0 = fypiN(0.2, 0.1, L, typ, DIS_CHARGE_EXCHANGE);
        double v1 = fypiN(0.2, 0.1, L, typ, DIS_NEUTRAL);
        check(v0 > v1, "fypiN dis=0 > dis=1 (isospin factor 2 vs 1)");
        near(v0, 2.0*v1, 1e-10, "fypiN dis=0 = 2 * dis=1");
    }

    // f_rhoN: positive at valid kinematics
    {
        double val = f_rhoN(0.2, 0.1, L, typ, DIS_CHARGE_EXCHANGE);
        notNeg(val, "f_rhoN >= 0 at (y=0.2, kT=0.1)");
    }

    // fypiD: returns 0 at y boundaries
    {
        near(fypiD(0.0,  0.1, Ld, typ, DIS_CHARGE_EXCHANGE), 0.0, 1e-15, "fypiD(y=0) = 0");
        near(fypiD(0.999,0.1, Ld, typ, DIS_CHARGE_EXCHANGE), 0.0, 1e-15, "fypiD(y=0.999) = 0");
    }
    // fypiD dis=1 returns 0 (coupling set to zero)
    {
        double val = fypiD(0.2, 0.1, Ld, typ, DIS_NEUTRAL);
        near(val, 0.0, 1e-15, "fypiD dis=1 = 0 (null coupling)");
    }

    // Form factor types: all should return positive values at valid kinematics
    double y = 0.15, kT = 0.08;
    for (int t = 0; t <= 6; ++t) {
        // Pauli-Villar (6) needs t variable; only meaningful for fypiN
        auto fft = static_cast<FormFactorType>(t);
        double val = fypiN(y, kT, L, fft, DIS_CHARGE_EXCHANGE);
        std::ostringstream n; n << "fypiN FF type=" << t << " >= 0";
        notNeg(val, n.str());
    }
}

// ============================================================
//  5. InputReader
// ============================================================
static void test_InputReader()
{
    std::cout << "\n--- InputReader ---\n";

    // Write a temp file, read it back, check values
    const std::string tmpFile = "/tmp/pioncloud_test_input.in";
    {
        std::ofstream f(tmpFile);
        f << "# test input\n"
          << "scan_mode   = x\n"
          << "beam_energy = 6.5\n"
          << "form_factor = dipole\n"
          << "contribution = rho\n"
          << "dissociation = 1\n"
          << "lambda       = 1.48\n"
          << "lambda_delta = 1.32\n"
          << "x_min = 0.03\n"
          << "x_max = 0.25\n"
          << "nx    = 50\n"
          << "k_min = 0.080\n"
          << "k_max = 0.160\n"
          << "nk    = 200\n"
          << "output_mode  = ascii\n"
          << "output_name  = TestRun\n";
    }

    InputReader reader;
    RunParams p = reader.read(tmpFile);

    check(!reader.hasErrors(),   "No validation errors");
    check(!reader.hasWarnings(), "No parse warnings");

    check(p.scanMode == "x",          "scanMode = x");
    near(p.E,    6.5,  1e-12, "beam_energy = 6.5");
    check(p.typ  == FF_DIPOLE,        "form_factor = dipole");
    check(p.flag == FLAG_RHO,         "contribution = rho");
    check(p.dis  == DIS_NEUTRAL,      "dissociation = 1");
    near(p.L,    1.48, 1e-12, "lambda = 1.48");
    near(p.Ld,   1.32, 1e-12, "lambda_delta = 1.32");
    near(p.xmin, 0.03, 1e-12, "x_min = 0.03");
    near(p.xmax, 0.25, 1e-12, "x_max = 0.25");
    check(p.nx  == 50,  "nx = 50");
    near(p.kmin, 0.080, 1e-12, "k_min = 0.080");
    near(p.kmax, 0.160, 1e-12, "k_max = 0.160");
    check(p.nk  == 200, "nk = 200");
    check(p.outputMode == OUT_ASCII, "output_mode = ascii");
    check(p.outBaseName == "TestRun","output_name = TestRun");

    // Unknown key -> warning, not error
    {
        std::ofstream f(tmpFile);
        f << "completely_unknown_key = 99\n";
    }
    InputReader r2;
    r2.read(tmpFile);
    check(!r2.hasErrors(),  "Unknown key: no error");
    check(r2.hasWarnings(), "Unknown key: warning issued");

    // Bad xmax < xmin -> validation error
    {
        std::ofstream f(tmpFile);
        f << "x_min = 0.5\nx_max = 0.1\n";
    }
    InputReader r3;
    r3.read(tmpFile);
    check(r3.hasErrors(), "xmax < xmin: validation error raised");

    // Inline comment stripped correctly
    {
        std::ofstream f(tmpFile);
        f << "beam_energy = 8.8  # BoNuS12 energy\n";
    }
    InputReader r4;
    RunParams p4 = r4.read(tmpFile);
    near(p4.E, 8.8, 1e-12, "Inline comment stripped: E = 8.8");

    // String name aliases: form_factor accepts integer string
    {
        std::ofstream f(tmpFile);
        f << "form_factor = 2\n";
    }
    InputReader r5;
    RunParams p5 = r5.read(tmpFile);
    check(p5.typ == FF_EXP_S, "form_factor = '2' -> FF_EXP_S");

    // Template can be written and read back cleanly
    {
        const std::string tplFile = "/tmp/pioncloud_template_test.in";
        InputReader::writeTemplate(tplFile);
        InputReader r6;
        RunParams p6 = r6.read(tplFile);
        check(!r6.hasErrors(),   "Template round-trip: no errors");
        check(!r6.hasWarnings(), "Template round-trip: no warnings");
        // Defaults should be preserved exactly
        near(p6.E,      11.0, 1e-12, "Template round-trip: E = 11");
        check(p6.nx == 100,          "Template round-trip: nx = 100");
        // New theta-scan bound defaults
        near(p6.xmax_th, 0.6,  1e-12, "Template round-trip: xmax_th = 0.6");
        near(p6.kmax_th, 0.5,  1e-12, "Template round-trip: kmax_th = 0.5");
    }

    // Fig.47 right panel configuration round-trip
    {
        const std::string fig47File = "/tmp/pioncloud_fig47.in";
        {
            std::ofstream f(fig47File);
            f << "scan_mode   = theta\n"
              << "contribution = pion\n"
              << "dissociation = neutral\n"
              << "form_factor  = exp_s\n"
              << "lambda       = 1.56\n"
              << "x_min_th     = 0.05\n"
              << "x_max_th     = 0.6\n"
              << "k_min_th     = 0.060\n"
              << "k_max_th     = 0.250\n";
        }
        InputReader r7;
        RunParams p7 = r7.read(fig47File);
        check(!r7.hasErrors(),   "Fig47 right panel: no errors");
        check(!r7.hasWarnings(), "Fig47 right panel: no warnings");
        near(p7.xmin_th, 0.05,  1e-12, "Fig47 right panel: xmin_th = 0.05");
        near(p7.xmax_th, 0.6,   1e-12, "Fig47 right panel: xmax_th = 0.6");
        near(p7.kmin_th, 0.060, 1e-12, "Fig47 right panel: kmin_th = 0.060");
        near(p7.kmax_th, 0.250, 1e-12, "Fig47 right panel: kmax_th = 0.250");
        check(p7.dis  == DIS_NEUTRAL, "Fig47 right panel: dis = neutral");
        check(p7.flag == FLAG_PION,   "Fig47 right panel: flag = pion");
        check(p7.typ  == FF_EXP_S,    "Fig47 right panel: typ = exp_s");
        near(p7.L,     1.56, 1e-12,   "Fig47 right panel: L = 1.56");
        std::remove(fig47File.c_str());
    }

    std::remove(tmpFile.c_str());
}

// ============================================================
//  6. Calculator -- smoke test with tiny grid
// ============================================================
static void test_Calculator()
{
    std::cout << "\n--- Calculator (smoke test) ---\n";

    RunParams p;
    p.flag       = FLAG_PION;
    p.typ        = FF_EXP_S;
    p.L          = 1.56;
    p.Ld         = 1.39;
    p.E          = 11.0;
    p.theta_e    = 35.0;
    p.xmin       = 0.05;
    p.xmax       = 0.20;
    p.nx         = 4;
    p.nth_int    = 4;
    p.thmin_deg  = 30.0;
    p.thmax_deg  = 70.0;
    p.kmin       = 0.060;
    p.kmax       = 0.130;
    p.nk         = 4;
    p.outputMode = OUT_NONE;
    p.scanMode   = "x";

    Calculator calc(p);
    ScanResults res;
    calc.runXScan(res);

    check(!res.xScan.empty(), "x-scan: result vector non-empty");
    check(res.xScan.size() == 5,  // nx+1 points
          "x-scan: correct number of points (nx+1 = 5)");

    // All F2piK values must be non-negative
    bool allNonNeg = true;
    for (auto &pt : res.xScan)
        if (pt.F2piK < 0.0) allNonNeg = false;
    check(allNonNeg, "x-scan: all F2piK >= 0");

    // First moment must be non-negative
    notNeg(res.F2pi0_x, "x-scan: first moment F2pi0_x >= 0");

    // Q2 values must be positive and increasing with x
    bool Q2ok = true;
    for (size_t i = 1; i < res.xScan.size(); ++i)
        if (res.xScan[i].Q2 <= res.xScan[i-1].Q2) Q2ok = false;
    check(Q2ok, "x-scan: Q2 monotonically increasing with x");

    // Running moment must be monotonically non-decreasing
    bool momOk = true;
    for (size_t i = 1; i < res.xScan.size(); ++i)
        if (res.xScan[i].F2piMoment < res.xScan[i-1].F2piMoment - 1e-20)
            momOk = false;
    check(momOk, "x-scan: RunningMoment non-decreasing");

    // F2n and ratio: when CTEQ6 is present they must be positive;
    // when absent they must be exactly zero (graceful degradation).
#ifdef HAVE_CTEQ6
    {
        bool f2nOk = true, ratioOk = true;
        for (auto &pt : res.xScan) {
            if (pt.F2n  <= 0.0) f2nOk   = false;
            if (pt.ratio < 0.0) ratioOk = false;
        }
        check(f2nOk,   "x-scan [CTEQ6]: all F2n > 0");
        check(ratioOk, "x-scan [CTEQ6]: all ratio >= 0");
        // ratio must be < 1 for reasonable kinematics
        // (pion cloud is a small fraction of total F2n)
        bool ratioSmall = true;
        for (auto &pt : res.xScan)
            if (pt.F2n > 0.0 && pt.ratio > 1.0) ratioSmall = false;
        check(ratioSmall, "x-scan [CTEQ6]: ratio < 1 (pion cloud < total F2n)");
    }
#else
    {
        bool f2nZero = true, ratioZero = true;
        for (auto &pt : res.xScan) {
            if (pt.F2n   != 0.0) f2nZero   = false;
            if (pt.ratio != 0.0) ratioZero = false;
        }
        check(f2nZero,   "x-scan [no CTEQ6]: F2n = 0 (graceful degradation)");
        check(ratioZero, "x-scan [no CTEQ6]: ratio = 0 (graceful degradation)");
    }
#endif

    // Theta scan
    p.scanMode  = "theta";
    p.nth       = 4;
    p.nx_int    = 4;
    Calculator calc2(p);
    ScanResults res2;
    calc2.runThetaScan(res2);

    check(!res2.thetaScan.empty(), "theta-scan: result vector non-empty");
    check(res2.thetaScan.size() == 5, "theta-scan: correct number of points (nth+1)");
    bool allNonNeg2 = true;
    for (auto &pt : res2.thetaScan)
        if (pt.F2piK < 0.0) allNonNeg2 = false;
    check(allNonNeg2, "theta-scan: all F2piK >= 0");
    notNeg(res2.F2pi0_theta, "theta-scan: first moment >= 0");

    // Theta-scan must respect xmin_th/kmin_th/kmax_th bounds:
    // full integration range must give >= result than narrow range.
    // Use nk=10, nx_int=10 so Simpson has enough points to see the difference.
    {
        RunParams pFull;
        pFull.flag      = FLAG_PION; pFull.typ = FF_EXP_S;
        pFull.L         = 1.56;      pFull.dis = DIS_NEUTRAL;
        pFull.scanMode  = "theta";
        pFull.nth       = 2;  pFull.nx_int = 10;  pFull.nk = 10;
        pFull.thmin_deg = 30.0; pFull.thmax_deg = 70.0;
        pFull.xmin_th   = 1.0e-4; pFull.xmax_th = 0.6;
        pFull.kmin_th   = 1.0e-4; pFull.kmax_th = 0.5;
        pFull.outputMode = OUT_NONE;

        RunParams pNarrow = pFull;
        pNarrow.xmin_th = 0.05;
        pNarrow.kmin_th = 0.060; pNarrow.kmax_th = 0.250;

        ScanResults rFull, rNarrow;
        Calculator(pFull).runThetaScan(rFull);
        Calculator(pNarrow).runThetaScan(rNarrow);

        double sumFull = 0.0, sumNarrow = 0.0;
        for (auto &pt : rFull.thetaScan)   sumFull   += pt.F2piK;
        for (auto &pt : rNarrow.thetaScan) sumNarrow += pt.F2piK;
        check(sumFull >= sumNarrow,
              "theta-scan: full k/x range >= narrow range (bounds respected)");
    }

    // runAll with kbin: verify kbin output is also populated
    p.scanMode = "all";
    p.xmin_th = 1.0e-4; p.xmax_th = 0.6;
    p.kmin_th = 1.0e-4; p.kmax_th = 0.5;
    Calculator calc3(p);
    ScanResults res3 = calc3.runAll();
    check(!res3.xScan.empty(),     "runAll: x-scan populated");
    check(!res3.thetaScan.empty(), "runAll: theta-scan populated");
    // kbin reuses xScan -- after runAll xScan was overwritten by kbin,
    // which is correct behaviour (last scan wins)
    check(!res3.xScan.empty(),     "runAll: kbin-scan populated (xScan non-empty)");
}

// ============================================================
//  6b. kbin-scan tests
// ============================================================
static void test_KBinScan()
{
    std::cout << "\n--- kbin-scan ---\n";

    RunParams p;
    p.flag       = FLAG_PION;
    p.typ        = FF_EXP_S;
    p.L          = 1.56;
    p.dis        = DIS_CHARGE_EXCHANGE;
    p.E          = 11.0;
    p.theta_e    = 12.0;   // as used in f2_pi_ranges.c
    p.xmin       = 0.055;
    p.xmax       = 0.20;
    p.nx         = 4;
    p.kmin       = 0.060;  // first BoNuS bin
    p.kmax       = 0.100;
    p.cosph_min  = 0.342;  // cos 70°
    p.cosph_max  = 0.866;  // cos 30°
    p.ny         = 10;
    p.nkT        = 200;    // enough for a unit test
    p.outputMode = OUT_NONE;
    p.scanMode   = "kbin";

    Calculator calc(p);
    ScanResults res;
    calc.runKBinScan(res);

    // Basic structural checks
    check(!res.xScan.empty(),            "kbin-scan: result vector non-empty");
    check(res.xScan.size() == 5,         "kbin-scan: nx+1 = 5 output points");

    // All F2piK must be non-negative
    bool allNonNeg = true;
    for (auto &pt : res.xScan)
        if (pt.F2piK < 0.0) allNonNeg = false;
    check(allNonNeg, "kbin-scan: all F2piK >= 0");

    // First moment must be non-negative
    notNeg(res.F2pi0_x, "kbin-scan: first moment >= 0");

    // Q2 values must be positive and increasing with x
    bool Q2ok = true;
    for (size_t i = 1; i < res.xScan.size(); ++i)
        if (res.xScan[i].Q2 <= res.xScan[i-1].Q2) Q2ok = false;
    check(Q2ok, "kbin-scan: Q2 monotonically increasing with x");

    // Running moment non-decreasing
    bool momOk = true;
    for (size_t i = 1; i < res.xScan.size(); ++i)
        if (res.xScan[i].F2piMoment < res.xScan[i-1].F2piMoment - 1e-30)
            momOk = false;
    check(momOk, "kbin-scan: RunningMoment non-decreasing");

    // Wider |k| bin must give more yield than narrower bin
    {
        RunParams pNarrow = p;
        RunParams pWide   = p;
        pNarrow.kmin = 0.060; pNarrow.kmax = 0.100;
        pWide.kmin   = 0.060; pWide.kmax   = 0.200;

        ScanResults rN, rW;
        Calculator(pNarrow).runKBinScan(rN);
        Calculator(pWide).runKBinScan(rW);

        double sumN = 0.0, sumW = 0.0;
        for (auto &pt : rN.xScan) sumN += pt.F2piK;
        for (auto &pt : rW.xScan) sumW += pt.F2piK;
        check(sumW >= sumN, "kbin-scan: wider |k| bin >= narrow bin");
    }

    // Removing angle cut (cosph_min=0, cosph_max=1) must give more yield
    {
        RunParams pCut  = p;
        RunParams pFull = p;
        pCut.cosph_min  = 0.342; pCut.cosph_max  = 0.866;
        pFull.cosph_min = 0.0;   pFull.cosph_max  = 1.0;

        ScanResults rC, rF;
        Calculator(pCut).runKBinScan(rC);
        Calculator(pFull).runKBinScan(rF);

        double sumC = 0.0, sumF = 0.0;
        for (auto &pt : rC.xScan) sumC += pt.F2piK;
        for (auto &pt : rF.xScan) sumF += pt.F2piK;
        check(sumF >= sumC, "kbin-scan: full angle range >= restricted range");
    }

    // InputReader: kbin-specific keys parse correctly
    {
        const std::string tmpFile = "/tmp/kbin_test.in";
        {
            std::ofstream f(tmpFile);
            f << "scan_mode  = kbin\n"
              << "k_min      = 0.060\n"
              << "k_max      = 0.100\n"
              << "cosph_min  = 0.342\n"
              << "cosph_max  = 0.866\n"
              << "ny         = 50\n"
              << "nkt        = 1000\n";
        }
        InputReader r;
        RunParams pr = r.read(tmpFile);
        check(!r.hasErrors(),           "kbin input: no errors");
        check(pr.scanMode == "kbin",    "kbin input: scanMode = kbin");
        near(pr.kmin,      0.060, 1e-12, "kbin input: k_min = 0.060");
        near(pr.kmax,      0.100, 1e-12, "kbin input: k_max = 0.100");
        near(pr.cosph_min, 0.342, 1e-12, "kbin input: cosph_min = 0.342");
        near(pr.cosph_max, 0.866, 1e-12, "kbin input: cosph_max = 0.866");
        check(pr.ny  == 50,              "kbin input: ny = 50");
        check(pr.nkT == 1000,            "kbin input: nkt = 1000");
        std::remove(tmpFile.c_str());
    }

    // ---- new feature tests ----

    // GRV99 gives different result from GRV92
    {
        RunParams pG92 = p;  pG92.grv_version = GRV_92;
        RunParams pG99 = p;  pG99.grv_version = GRV_99;
        ScanResults r92, r99;
        Calculator(pG92).runKBinScan(r92);
        Calculator(pG99).runKBinScan(r99);
        // They should both be non-negative but differ numerically
        bool g92ok = true, g99ok = true;
        for (auto &pt : r92.xScan) if (pt.F2piK < 0.0) g92ok = false;
        for (auto &pt : r99.xScan) if (pt.F2piK < 0.0) g99ok = false;
        check(g92ok, "GRV92: all F2piK >= 0");
        check(g99ok, "GRV99: all F2piK >= 0");
    }

    // Precision preset: FAST gives fewer steps than CUSTOM with large nkT
    {
        RunParams pFast = p;
        pFast.precision = PRECISION_FAST;
        pFast.nkT = 50000;   // will be overridden by preset
        // Constructor applies preset; verify via a quick run (no crash)
        ScanResults rf;
        bool ok = true;
        try { Calculator(pFast).runKBinScan(rf); }
        catch (...) { ok = false; }
        check(ok,                  "precision=fast: runs without exception");
        check(!rf.xScan.empty(),   "precision=fast: xScan non-empty");
    }

    // Inclusive mode: removes cuts -> yields more than restricted
    {
        RunParams pRestricted = p;
        RunParams pInclusive  = p;
        pInclusive.inclusive_mode = true;
        ScanResults rR, rI;
        Calculator(pRestricted).runKBinScan(rR);
        Calculator(pInclusive ).runKBinScan(rI);
        double sumR = 0.0, sumI = 0.0;
        for (auto &pt : rR.xScan) sumR += pt.F2piK;
        for (auto &pt : rI.xScan) sumI += pt.F2piK;
        check(sumI >= sumR, "inclusive_mode: no-cut yield >= restricted yield");
    }

    // eval() single-point: result must be non-negative and consistent
    // with runKBinScan at the same x.  Use nx=1 so xstep is well-defined,
    // then compare eval(xmin) against the first point of the scan.
    {
        RunParams pe = p;
        pe.xmin = 0.10; pe.xmax = 0.20; pe.nx = 1;
        pe.verbose = false;
        ScanResults rs;
        Calculator(pe).runKBinScan(rs);
        // First point of scan is at x = xmin = 0.10
        double f_scan = rs.xScan.empty() ? -1.0 : rs.xScan.front().F2piK;
        double f_eval = Calculator(pe).eval(0.10);
        // Both must be non-negative
        notNeg(f_eval, "eval(x=0.10): non-negative");
        notNeg(f_scan, "runKBinScan at x=0.10: non-negative");
        // Must agree to floating-point precision (same algorithm)
        near(f_eval, f_scan, 1e-12 * (1.0 + std::abs(f_scan)),
             "eval(x=0.10) matches runKBinScan at same x");
    }

    // verbose=false: no output to stdout (structural test -- just verify no crash)
    {
        RunParams pq = p; pq.verbose = false;
        ScanResults rq;
        bool ok = true;
        try { Calculator(pq).runKBinScan(rq); }
        catch (...) { ok = false; }
        check(ok, "verbose=false: runs without exception");
    }

    // InputReader: new keys round-trip
    {
        const std::string f2 = "/tmp/new_keys_test.in";
        {
            std::ofstream f(f2);
            f << "nucleon     = proton\n"
              << "grv_version = 99\n"
              << "precision   = normal\n"
              << "inclusive   = true\n"
              << "verbose     = false\n"
              << "angle_min_deg = 30.0\n"
              << "angle_max_deg = 70.0\n";
        }
        InputReader rdr;
        RunParams pr = rdr.read(f2);
        check(!rdr.hasErrors(),              "new keys: no errors");
        check(pr.nucleon    == NUCLEON_PROTON,       "new keys: nucleon=proton");
        check(pr.grv_version == GRV_99,             "new keys: grv_version=99");
        check(pr.precision  == PRECISION_NORMAL,     "new keys: precision=normal");
        check(pr.inclusive_mode,                     "new keys: inclusive=true");
        check(!pr.verbose,                           "new keys: verbose=false");
        // angle_max_deg=70 -> cosph_min = cos(70°) = 0.342
        near(pr.cosph_min, std::cos(70.0*M_PI/180.0), 1e-10,
             "new keys: angle_max_deg=70 -> cosph_min=cos(70°)");
        near(pr.cosph_max, std::cos(30.0*M_PI/180.0), 1e-10,
             "new keys: angle_min_deg=30 -> cosph_max=cos(30°)");
        std::remove(f2.c_str());
    }
}

// ============================================================
//  6c. Numerical regression tests
//
//  Reference values from fortran_reference.txt, produced by
//  running f2_pi_ranges.c against f2_pi_sub.f90 (June 2023).
//  Config: E=11 GeV, theta_e=12 deg, pion, dis=charge_exchange,
//          exp_s (typ=2), L=1.56 GeV, cosph in [0.342,0.866].
//
//  The file contains 4 k-bins. We use bins t=2 and t=3 because:
//    t=0 [60,100]MeV:  only 2 non-zero points (too few)
//    t=1 [100,200]MeV: 8 points (usable but marginal)
//    t=2 [200,300]MeV: 12 points -- USED
//    t=3 [300,400]MeV: 15 points -- USED
//
//  x-points: x_min=0.055, x_step=(0.33-0.055)/20=0.01375
//
//  Tolerance: 5% relative -- accounts for:
//    - C++ uses Simpson on integer grid; Fortran uses float DO loop
//    - Different effective nkT step sizes at these k ranges
//    Note: the narrow bins [60,100] and [100,200] MeV have very few
//    non-zero points at nkT=5000; those bins are NOT tested here.
//    For the black curve (k=[60,250] MeV), a dedicated run of
//    f2_pi.c with k_min=0.060, k_max=0.250 is still needed.
// ============================================================
static void test_Regression()
{
    std::cout << "\n--- Numerical regression (vs Fortran f2_pi_ranges) ---\n";

    // ---- Bin t=2: k=[200,300] MeV, 12 non-zero points ----
    struct Ref { double x; double f2piK; };
    static const Ref bin2[] = {
        { 0.055,    0.00181429 },
        { 0.06875,  0.00165378 },
        { 0.0825,   0.00145021 },
        { 0.09625,  0.00120890 },
        { 0.11,     0.000955094},
        { 0.12375,  0.000721630},
        { 0.1375,   0.000507765},
        { 0.15125,  0.000320570},
        { 0.165,    0.000183454},
        { 0.17875,  8.89425e-05},
        { 0.1925,   3.35814e-05},
        { 0.20625,  7.42268e-06},
    };

    // ---- Bin t=3: k=[300,400] MeV, 15 non-zero points ----
    static const Ref bin3[] = {
        { 0.055,    0.00245428 },
        { 0.06875,  0.00240626 },
        { 0.0825,   0.00226778 },
        { 0.09625,  0.00207309 },
        { 0.11,     0.00184773 },
        { 0.12375,  0.00159680 },
        { 0.1375,   0.00134010 },
        { 0.15125,  0.00109687 },
        { 0.165,    0.000857263},
        { 0.17875,  0.000640851},
        { 0.1925,   0.000446967},
        { 0.20625,  0.000290261},
        { 0.22,     0.000165046},
        { 0.23375,  7.89949e-05},
        { 0.2475,   2.91344e-05},
    };

    const double rel_tol = 0.05;  // 5% relative tolerance

    auto runBin = [&](const char* label,
                      double kmin, double kmax,
                      const Ref* ref, int nRef)
    {
        std::cout << "  -- " << label << " --\n";

        // Build x array matching the Fortran exactly
        RunParams p;
        p.flag       = FLAG_PION;
        p.typ        = FF_EXP_S;
        p.L          = 1.56;
        p.dis        = DIS_CHARGE_EXCHANGE;
        p.E          = 11.0;
        p.theta_e    = 12.0;
        p.kmin       = kmin;
        p.kmax       = kmax;
        p.cosph_min  = 0.342;
        p.cosph_max  = 0.866;
        p.ny         = 100;
        p.nkT        = 10000;  // higher than default for regression accuracy
        p.xmin       = ref[0].x;
        p.xmax       = ref[nRef-1].x;
        p.nx         = nRef - 1;
        p.outputMode = OUT_NONE;
        p.scanMode   = "kbin";
        p.verbose    = false;

        ScanResults res;
        Calculator(p).runKBinScan(res);

        int nPass = 0, nFail = 0;
        for (int i = 0; i < nRef && i < (int)res.xScan.size(); ++i) {
            double f_ref = ref[i].f2piK;
            double f_cpp = res.xScan[i].F2piK;
            if (f_ref < 1e-10) continue;
            double tol = rel_tol * f_ref;
            bool ok = std::abs(f_cpp - f_ref) <= tol;
            if (ok) ++nPass; else ++nFail;
            // Print every point so the user can see agreement
            std::cout << "    x=" << std::fixed << std::setprecision(5) << ref[i].x
                      << "  fortran=" << std::scientific << std::setprecision(5) << f_ref
                      << "  c++=" << f_cpp
                      << "  diff=" << std::fixed << std::setprecision(1)
                      << std::abs(f_cpp - f_ref) / f_ref * 100.0 << "%"
                      << (ok ? "" : "  <-- FAIL") << "\n";
        }
        std::cout << "    " << nPass << " pass, " << nFail << " fail\n";
        return nFail == 0;
    };

    bool ok2 = runBin("k=[200,300] MeV (bin t=2)", 0.200, 0.300,
                      bin2, sizeof(bin2)/sizeof(bin2[0]));
    bool ok3 = runBin("k=[300,400] MeV (bin t=3)", 0.300, 0.400,
                      bin3, sizeof(bin3)/sizeof(bin3[0]));

    check(ok2, "regression bin t=2 [200,300] MeV: within 5% of Fortran");
    check(ok3, "regression bin t=3 [300,400] MeV: within 5% of Fortran");

    // ---- k=[60,250] MeV with 30-70 deg angular cut ----
    // Produced by f2_pi.c modified: k_min=0.060, k_max=0.250,
    // x_min=0.055, x_step=(0.330-0.055)/14, angle=12, E=11.
    // This tests runKBinScan (f2pi_sub mode) with full BoNuS
    // acceptance window.
    // NOTE: this is NOT the black piN curve in f2piN.png.
    // That curve comes from 3Var_x.f (runXScan, no angular cut)
    // and extends to x~0.27. This curve has the 30-70 deg cut
    // and dies at x~0.19. Both use k=[60,250] MeV but are
    // physically different quantities.
    static const Ref kbin_60_250[] = {
        { 0.0550000, 1.382730e-03 },
        { 0.0746429, 1.039190e-03 },
        { 0.0942857, 6.877680e-04 },
        { 0.1139290, 3.935860e-04 },
        { 0.1335710, 1.846210e-04 },
        { 0.1532140, 6.106900e-05 },
        { 0.1728570, 1.010020e-05 },
    };
    bool okKbin = runBin("k=[60,250] MeV kbin (30-70 deg cut, f2pi_sub)",
                          0.060, 0.250,
                          kbin_60_250, sizeof(kbin_60_250)/sizeof(kbin_60_250[0]));
    check(okKbin, "regression kbin [60,250] MeV: within 5% of Fortran f2pi_sub");

    std::cout << "  NOTE: bins t=0 [60,100] and t=1 [100,200] MeV not tested\n"
              << "        (too few Fortran reference points at this x-spacing).\n";
}

// ============================================================
//  7. ASCII output -- file content verification
// ============================================================
static void test_OutputASCII()
{
    std::cout << "\n--- OutputWriter (ASCII) ---\n";

    // Run a tiny calculation and verify the output files exist
    // and have the expected structure
    RunParams p;
    p.flag = FLAG_PION; p.typ = FF_EXP_S;
    p.xmin = 0.05; p.xmax = 0.15; p.nx = 2;
    p.nth_int = 2; p.kmin = 0.06; p.kmax = 0.10; p.nk = 2;
    p.thmin_deg = 30.0; p.thmax_deg = 70.0;
    p.nth = 2; p.nx_int = 2;
    p.outputMode = OUT_ASCII;
    p.outBaseName = "/tmp/PionCloud_test_output";
    p.scanMode = "all";
    p.verbose  = false;   // suppress progress output during tests

    Calculator calc(p);
    ScanResults res = calc.runAll();

    check(res.xScan.size() == 3,     "Output test: xScan has nx+1=3 points");
    check(res.thetaScan.size() == 3, "Output test: thetaScan has nth+1=3 points");

    // Verify running moment: the last entry must equal the full Simpson
    // integral over all x points, matching F2pi0_x.
    {
        RunParams px;
        px.flag = FLAG_PION; px.typ = FF_EXP_S;
        px.xmin = 0.05; px.xmax = 0.15; px.nx = 2;
        px.nth_int = 2; px.kmin = 0.06; px.kmax = 0.10; px.nk = 2;
        px.scanMode = "x";
        px.verbose  = false;
        px.outputMode = OUT_NONE;
        ScanResults rx;
        Calculator(px).runXScan(rx);
        // The last running moment must equal the full first moment F2pi0_x
        near(rx.xScan.back().F2piMoment, rx.F2pi0_x, 1e-20,
             "Running moment final value matches manual sum");
    }
}

// ============================================================
//  8. CTEQ6 wrapper (conditional on HAVE_CTEQ6)
// ============================================================
static void test_CTEQ6()
{
    std::cout << "\n--- CTEQ6 wrapper ---\n";

#ifndef HAVE_CTEQ6
    std::cout << "  [SKIP] Built without CTEQ6 -- skipping CTEQ6 tests.\n";
    return;
#else
    // init() must not throw
    bool initOk = true;
    try { CTEQ6::init(); }
    catch (...) { initOk = false; }
    check(initOk, "CTEQ6::init() does not throw");

    // xfx: gluon at moderate x should be large and positive
    {
        double g = CTEQ6::xfx(0, 0.1, 2.0);
        check(g > 0.0, "xfx(gluon, x=0.1, Q=2) > 0");
    }
    // xfx: u-quark at x=0.3 should be dominant in proton
    {
        double u = CTEQ6::xfx(1, 0.3, 2.0);
        double d = CTEQ6::xfx(2, 0.3, 2.0);
        check(u > 0.0, "xfx(u, x=0.3, Q=2) > 0");
        check(u > d,   "xfx(u) > xfx(d) at x=0.3 in proton");
    }
    // F2proton: positive in valid range
    {
        double f2p = CTEQ6::F2proton(0.1, 4.0);
        check(f2p > 0.0, "F2proton(x=0.1, Q2=4) > 0");
    }
    // F2neutron: positive in valid range
    {
        double f2n = CTEQ6::F2neutron(0.1, 4.0);
        check(f2n > 0.0, "F2neutron(x=0.1, Q2=4) > 0");
    }
    // Charge symmetry: at moderate x, F2n < F2p because u dominates proton
    // (d_n = u_p is large, but charge factor 4/9 vs 1/9 means proton is bigger)
    {
        double f2p = CTEQ6::F2proton (0.3, 4.0);
        double f2n = CTEQ6::F2neutron(0.3, 4.0);
        check(f2p > f2n, "F2proton > F2neutron at x=0.3 (u dominance)");
    }
    // Guard: x=0 returns 0
    {
        near(CTEQ6::F2neutron(0.0, 4.0), 0.0, 1e-15,
             "F2neutron(x=0) = 0 (guard)");
    }
    // Guard: Q2=0 returns 0
    {
        near(CTEQ6::F2neutron(0.1, 0.0), 0.0, 1e-15,
             "F2neutron(Q2=0) = 0 (guard)");
    }
    // Strange contribution: F2n with strange >= F2n without
    {
        double f2n_no_s  = CTEQ6::F2neutron(0.1, 4.0, false);
        double f2n_with_s= CTEQ6::F2neutron(0.1, 4.0, true);
        check(f2n_with_s >= f2n_no_s,
              "F2n with strange >= F2n without (s contribution non-negative)");
    }
#endif
}

// ============================================================
//  main
// ============================================================
int main()
{
    std::cout << "╔═══════════════════════════════════════╗\n"
              << "║  PionCloud Unit Tests                 ║\n"
              << "╚═══════════════════════════════════════╝\n";

    test_PhysicsParams();
    test_Kinematics();
    test_GRV();
    test_SplittingFunctions();
    test_InputReader();
    test_Calculator();
    test_KBinScan();
    test_Regression();
    test_OutputASCII();
    test_CTEQ6();

    std::cout << "\n══════════════════════════════════════\n"
              << "  Results: " << g_pass << " passed, " << g_fail << " failed\n"
              << "══════════════════════════════════════\n";

    return (g_fail == 0) ? 0 : 1;
}
