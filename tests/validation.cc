// ============================================================
//  tests/validation.cc   (new in v1.4)
//  Regression test against the PUBLISHED TDIS curves.
//
//  Runs the configurations of docs/validation/make_runs.sh with the
//  Calculator and compares with the curves extracted from the vector
//  figures of PR12-15-006 (docs/validation/ref_data, see README there).
//
//  Criteria (on the "bulk": published curve within x100 of its maximum,
//  where the comparison is meaningful; near the kinematic end points both
//  curves drop steeply to zero):
//    tagged curves     : -1% < median deviation < +3%,  max |dev| < 8%
//    total F2^(piN)    : |median deviation| < 0.5%  (total_charge_states)
//  Reference values (v1.4, grid method): medians +0.2..+1.7%, max <= 6.9%.
//  The published curves sit ~+1% above PionCloud because the original
//  Fortran Simpson loops give the last grid point weight 2 instead of 1.
//
//  Output and pass criterion follow tests/test_all.cc ("Results: N passed,
//  0 failed" is required by ctest).
// ============================================================
#include "Calculator.hh"
#include "PhysicsParams.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace PionCloud;

#ifndef PIONCLOUD_REF_DIR
#error "PIONCLOUD_REF_DIR must point to docs/validation/ref_data"
#endif

static int g_pass = 0, g_fail = 0;

static std::vector<std::pair<double,double>> readRef(const std::string &name)
{
    std::ifstream f(std::string(PIONCLOUD_REF_DIR) + "/" + name + ".dat");
    std::vector<std::pair<double,double>> v;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line); double a, b;
        if (is >> a >> b) v.emplace_back(a, b);
    }
    std::sort(v.begin(), v.end());
    return v;
}

// log-interpolate (xs, ys) at x0 (ys > 0)
static double logInterp(const std::vector<double> &xs, const std::vector<double> &ys, double x0)
{
    for (size_t i = 1; i < xs.size(); ++i)
        if (x0 <= xs[i]) {
            const double t = (x0 - xs[i-1]) / (xs[i] - xs[i-1]);
            return std::exp((1-t)*std::log(ys[i-1]) + t*std::log(ys[i]));
        }
    return ys.back();
}

static void compare(const std::string &label, const std::string &ref,
                    const std::vector<double> &xs_in, const std::vector<double> &ys_in,
                    double medLo, double medHi, double maxDev)
{
    std::vector<double> xs, ys;
    for (size_t i = 0; i < xs_in.size(); ++i)
        if (ys_in[i] > 0.0) { xs.push_back(xs_in[i]); ys.push_back(ys_in[i]); }
    auto R = readRef(ref);
    double rmax = 0.0; for (auto &p : R) rmax = std::max(rmax, p.second);
    std::vector<double> dev;
    for (auto &p : R)
        if (p.first >= xs.front() && p.first <= xs.back() && p.second > 1e-2 * rmax)
            dev.push_back(logInterp(xs, ys, p.first) / p.second - 1.0);
    bool ok = !R.empty() && dev.size() >= 3;
    double med = 0.0, worst = 0.0;
    if (ok) {
        std::vector<double> d = dev; std::sort(d.begin(), d.end());
        med = d[d.size()/2];
        for (double e : dev) worst = std::max(worst, std::abs(e));
        ok = (med > medLo && med < medHi && worst < maxDev);
    }
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << std::left << std::setw(30) << label
              << std::right << std::fixed << std::setprecision(2)
              << "  median " << std::setw(6) << 100*med << "%  max|dev| " << std::setw(5) << 100*worst << "%  ("
              << dev.size() << " pts)\n";
    (ok ? g_pass : g_fail)++;
}

static RunParams base()
{
    RunParams p;                          // HSS exp_s defaults (v1.4)
    p.E = 11.0; p.theta_e = 35.0;
    p.typ = FF_EXP_S; p.L = 1.56; p.Ld = 1.39;
    p.outputMode = OUT_NONE; p.verbose = false;
    return p;
}

int main()
{
    std::cout << "PionCloud validation against PR12-15-006 (Figs. 5, 6, 37, 46, 47)\n";
    const double TLo = -0.01, THi = 0.03, TMax = 0.08;    // tagged curves
    const double ILo = -0.005, IHi = 0.005, IMax = 0.01;  // totals

    // ---- Fig. 47: theta scans ----
    for (int panel = 0; panel < 2; ++panel)
        for (ContributionFlag fl : {FLAG_PION, FLAG_RHO}) {
            RunParams p = base();
            p.flag = fl; p.dis = DIS_NEUTRAL;
            p.thmin_deg = 0; p.thmax_deg = 90; p.nth = 90; p.nx_int = 200; p.nk = 200;
            if (panel == 1) { p.xmin_th = 0.05; p.kmin_th = 0.060; p.kmax_th = 0.250; }
            ScanResults r; Calculator(p).runThetaScan(r);
            std::vector<double> xs, ys;
            for (auto &pt : r.thetaScan) { xs.push_back(pt.theta_deg); ys.push_back(pt.F2piK); }
            const std::string side = panel ? "right" : "left";
            compare("Fig.47 " + side + (fl == FLAG_PION ? " pi0 p" : " rho0 p"),
                    "fig47_" + side + (fl == FLAG_PION ? "_red" : "_blue_dash"),
                    xs, ys, TLo, THi, TMax);
        }

    // ---- Fig. 46: x scan, all angles, Delta k = [0,500] MeV ----
    for (ContributionFlag fl : {FLAG_PION, FLAG_RHO}) {
        RunParams p = base();
        p.flag = fl; p.dis = DIS_NEUTRAL;
        p.xmin = 0.005; p.xmax = 0.4; p.nx = 80;
        p.thmin_deg = 0; p.thmax_deg = 180; p.nth_int = 200;
        p.kmin = 1e-4; p.kmax = 0.5; p.nk = 200;
        ScanResults r; Calculator(p).runXScan(r);
        std::vector<double> xs, ys;
        for (auto &pt : r.xScan) { xs.push_back(pt.x); ys.push_back(pt.F2piK); }
        compare(std::string("Fig.46 ") + (fl == FLAG_PION ? "pi0 p" : "rho0 p"),
                fl == FLAG_PION ? "fig46_red" : "fig46_blue_dash", xs, ys, TLo, THi, TMax);
    }

    // ---- kbin figures ----
    struct K { const char *label, *ref; ContributionFlag fl; DisChannel dis;
               double kmin, kmax, cmin, xmax; bool total; NucleonType nuc; };
    const K ks[] = {
        {"Fig.5 |k| 100-200", "fig5_main_red_dash",   FLAG_PION, DIS_NEUTRAL, 0.100, 0.200, 0.342, 0.30, false, NUCLEON_PROTON},
        {"Fig.5 |k| 200-300", "fig5_main_green_dash", FLAG_PION, DIS_NEUTRAL, 0.200, 0.300, 0.342, 0.30, false, NUCLEON_PROTON},
        {"Fig.5 |k| 300-400", "fig5_main_blue",       FLAG_PION, DIS_NEUTRAL, 0.300, 0.400, 0.342, 0.30, false, NUCLEON_PROTON},
        {"Fig.5 total F2(pi p)", "fig5_main_violet_dash", FLAG_PION, DIS_NEUTRAL, -1, -1, 0, 0.30, true, NUCLEON_PROTON},
        {"Fig.6 piN",        "fig6_main_black",      FLAG_PION,     DIS_CHARGE_EXCHANGE, 0.060, 0.250, 0.0, 0.20, false, NUCLEON_NEUTRON},
        {"Fig.6 rhoN",       "fig6_main_red_dash",   FLAG_RHO,      DIS_CHARGE_EXCHANGE, 0.060, 0.250, 0.0, 0.20, false, NUCLEON_NEUTRON},
        {"Fig.6 piDelta",    "fig6_main_green_dash", FLAG_PI_DELTA, DIS_CHARGE_EXCHANGE, 0.060, 0.250, 0.0, 0.20, false, NUCLEON_NEUTRON},
        {"Fig.6 total F2(piN)", "fig6_main_violet_dash", FLAG_PION, DIS_CHARGE_EXCHANGE, -1, -1, 0, 0.30, true, NUCLEON_NEUTRON},
        {"Fig.37 piN",       "fig37_main_black",     FLAG_PION,     DIS_CHARGE_EXCHANGE, 0.250, 0.400, 0.342, 0.30, false, NUCLEON_NEUTRON},
        {"Fig.37 rhoN",      "fig37_main_red_dash",  FLAG_RHO,      DIS_CHARGE_EXCHANGE, 0.250, 0.400, 0.342, 0.30, false, NUCLEON_NEUTRON},
        {"Fig.37 piDelta",   "fig37_main_green_dash",FLAG_PI_DELTA, DIS_CHARGE_EXCHANGE, 0.250, 0.400, 0.342, 0.30, false, NUCLEON_NEUTRON},
    };
    // v1.4: each kbin figure is checked with BOTH integration methods --
    // the original Fortran grid (with which the figures were made) and the
    // new default exact method.
    for (KBinMethod meth : {KBIN_GRID, KBIN_EXACT})
    for (const K &k : ks) {
        RunParams p = base();
        p.kbin_method = meth;
        p.flag = k.fl; p.dis = k.dis; p.nucleon = k.nuc;
        p.xmin = 0.05; p.xmax = k.xmax; p.nx = 50; p.ny = 100; p.nkT = 5000;
        if (k.total) { p.inclusive_mode = true; p.total_charge_states = true; }
        else { p.kmin = k.kmin; p.kmax = k.kmax; p.cosph_min = k.cmin; p.cosph_max = 0.866; }
        ScanResults r; Calculator(p).runKBinScan(r);
        std::vector<double> xs, ys;
        for (auto &pt : r.xScan) { xs.push_back(pt.x); ys.push_back(pt.F2piK); }
        const std::string lab = std::string(k.label) + (meth == KBIN_GRID ? " [grid]" : " [exact]");
        if (k.total) compare(lab, k.ref, xs, ys, ILo, IHi, IMax);
        else         compare(lab, k.ref, xs, ys, TLo, THi, TMax);
    }

    std::cout << "  Results: " << g_pass << " passed, " << g_fail << " failed\n";
    return g_fail == 0 ? 0 : 1;
}
