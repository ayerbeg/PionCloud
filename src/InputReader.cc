// ============================================================
//  InputReader.cc
//  Pion Cloud Model -- input file parser (implementation)
// ============================================================

#include "InputReader.hh"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PionCloud {

// ============================================================
//  String helpers
// ============================================================
std::string InputReader::trim(const std::string &s)
{
    auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string InputReader::toLower(const std::string &s)
{
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return r;
}

// ============================================================
//  read  -- open file, return populated RunParams
// ============================================================
RunParams InputReader::read(const std::string &filename)
{
    RunParams p;          // starts with all defaults
    readInto(filename, p);
    return p;
}

// ============================================================
//  readInto  -- parse file, merge into caller's RunParams
// ============================================================
void InputReader::readInto(const std::string &filename, RunParams &params)
{
    warnings_.clear();
    errors_.clear();

    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("InputReader: cannot open '" + filename + "'");

    std::string line;
    int lineNo = 0;
    while (std::getline(f, line)) {
        ++lineNo;
        parseLine(line, lineNo, params);
    }

    // Run validation after all keys are read
    errors_ = validate(params);
}

// ============================================================
//  parseLine  -- handle one line of the input file
//
//  Accepted formats:
//    key = value        # inline comment
//    # full-line comment
//    [section]          (informational only, ignored)
//    blank line
// ============================================================
void InputReader::parseLine(const std::string &rawLine,
                             int lineNo, RunParams &p)
{
    std::string line = trim(rawLine);

    // skip blank and comment lines
    if (line.empty() || line[0] == '#' || line[0] == ';') return;
    // skip section headers  [...]
    if (line[0] == '[') return;

    // strip inline comment
    auto cpos = line.find('#');
    if (cpos != std::string::npos) line = trim(line.substr(0, cpos));

    // split on '='
    auto eq = line.find('=');
    if (eq == std::string::npos) {
        warnings_.push_back({lineNo, "No '=' found, line ignored: '" + line + "'"});
        return;
    }

    std::string key = toLower(trim(line.substr(0, eq)));
    std::string val =         trim(line.substr(eq + 1));

    if (key.empty() || val.empty()) {
        warnings_.push_back({lineNo, "Empty key or value, line ignored: '" + line + "'"});
        return;
    }

    // ---- helper lambdas ----
    auto asDouble = [&](double &dst) {
        try { dst = std::stod(val); }
        catch (...) {
            warnings_.push_back({lineNo,
                "Cannot parse '" + val + "' as double for key '" + key + "'"});
        }
    };
    auto asInt = [&](int &dst) {
        try { dst = std::stoi(val); }
        catch (...) {
            warnings_.push_back({lineNo,
                "Cannot parse '" + val + "' as int for key '" + key + "'"});
        }
    };

    // ---- key dispatch ----
    // [Beam Kinematics]
    if      (key == "e" || key == "beam_energy")        asDouble(p.E);
    else if (key == "theta_e" || key == "electron_angle") asDouble(p.theta_e);

    // [Physics - Cutoffs]
    else if (key == "lambda" || key == "l")             asDouble(p.L);
    else if (key == "lambda_delta" || key == "ld")      asDouble(p.Ld);

    // [Physics - Form Factor]
    else if (key == "form_factor" || key == "typ") {
        // Accept either integer code or name string
        std::string vl = toLower(val);
        if      (vl == "0" || vl == "monopole")          p.typ = FF_MONOPOLE;
        else if (vl == "1" || vl == "dipole")             p.typ = FF_DIPOLE;
        else if (vl == "2" || vl == "exp_s" ||
                 vl == "exponential_s")                   p.typ = FF_EXP_S;
        else if (vl == "3" || vl == "cov_dipole" ||
                 vl == "covariant_dipole")                p.typ = FF_COV_DIPOLE;
        else if (vl == "4" || vl == "dipole_sch" ||
                 vl == "dipole_s_channel")                p.typ = FF_DIPOLE_S_CH;
        else if (vl == "5" || vl == "exp_t" ||
                 vl == "exponential_t")                   p.typ = FF_EXP_T;
        else if (vl == "6" || vl == "pauli_villar" ||
                 vl == "paulivillar" || vl == "pv")       p.typ = FF_PAULI_VILLAR;
        else warnings_.push_back({lineNo,
            "Unknown form_factor value '" + val + "', keeping current setting"});
    }

    // [Physics - Contribution flag]
    else if (key == "contribution" || key == "flag") {
        std::string vl = toLower(val);
        if      (vl == "0" || vl == "pion"  || vl == "pi_n")       p.flag = FLAG_PION;
        else if (vl == "1" || vl == "rho"   || vl == "rho_n")      p.flag = FLAG_RHO;
        else if (vl == "2" || vl == "pi_delta" || vl == "pidelta")  p.flag = FLAG_PI_DELTA;
        else warnings_.push_back({lineNo,
            "Unknown contribution value '" + val + "', keeping current setting"});
    }

    // [Physics - Dissociation channel]
    else if (key == "dissociation" || key == "dis") {
        std::string vl = toLower(val);
        if      (vl == "0" || vl == "charge_exchange" ||
                 vl == "charge")                          p.dis = DIS_CHARGE_EXCHANGE;
        else if (vl == "1" || vl == "neutral")            p.dis = DIS_NEUTRAL;
        else warnings_.push_back({lineNo,
            "Unknown dissociation value '" + val + "', keeping current setting"});
    }

    // [x-Scan]
    else if (key == "x_min" || key == "xmin")   asDouble(p.xmin);
    else if (key == "x_max" || key == "xmax")   asDouble(p.xmax);
    else if (key == "nx" || key == "n_x")       asInt   (p.nx);

    // [theta-Scan]
    else if (key == "theta_min" || key == "thmin") asDouble(p.thmin_deg);
    else if (key == "theta_max" || key == "thmax") asDouble(p.thmax_deg);
    else if (key == "nth" || key == "n_theta")     asInt   (p.nth);
    // theta-scan internal x integration bounds
    // Left panel (default): x_min_th~0, x_max_th=0.6
    // Right panel (Fig.47): x_min_th=0.05, x_max_th=0.6
    else if (key == "x_min_th" || key == "xmin_th") asDouble(p.xmin_th);
    else if (key == "x_max_th" || key == "xmax_th") asDouble(p.xmax_th);
    // theta-scan internal |k| integration bounds
    // Left panel (default): k_min_th~0, k_max_th=0.5
    // Right panel (Fig.47): k_min_th=0.060, k_max_th=0.250
    else if (key == "k_min_th" || key == "kmin_th") asDouble(p.kmin_th);
    else if (key == "k_max_th" || key == "kmax_th") asDouble(p.kmax_th);

    // [Hadron kinematics - from subroutine]
    else if (key == "ph" || key == "hadron_momentum")   asDouble(p.pH);
    else if (key == "alpha1" || key == "hadron_angle_min") asDouble(p.alpha1_deg);
    else if (key == "alpha2" || key == "hadron_angle_max") asDouble(p.alpha2_deg);

    // [|k| integration]
    else if (key == "k_min" || key == "kmin")   asDouble(p.kmin);
    else if (key == "k_max" || key == "kmax")   asDouble(p.kmax);
    else if (key == "nk" || key == "n_k")       asInt   (p.nk);

    // [Output]
    else if (key == "output" || key == "output_mode") {
        std::string vl = toLower(val);
        if      (vl == "none")          p.outputMode = OUT_NONE;
        else if (vl == "ascii")         p.outputMode = OUT_ASCII;
        else if (vl == "root")          p.outputMode = OUT_ROOT;
        else if (vl == "all" ||
                 vl == "ascii+root" ||
                 vl == "root+ascii")    p.outputMode = OUT_ALL;
        else warnings_.push_back({lineNo,
            "Unknown output value '" + val + "', keeping current setting"});
    }
    else if (key == "output_name" || key == "basename" ||
             key == "out_name")         p.outBaseName = val;

    // [Physics - nucleon type]
    else if (key == "nucleon") {
        std::string vl = toLower(val);
        if      (vl == "neutron" || vl == "0") p.nucleon = NUCLEON_NEUTRON;
        else if (vl == "proton"  || vl == "1") p.nucleon = NUCLEON_PROTON;
        else warnings_.push_back({lineNo,
            "Unknown nucleon '" + val + "', expected neutron|proton"});
    }
    // [Physics - GRV pion PDF version]
    else if (key == "grv_version" || key == "grv") {
        std::string vl = toLower(val);
        if      (vl == "92" || vl == "grv92") p.grv_version = GRV_92;
        else if (vl == "99" || vl == "grv99") p.grv_version = GRV_99;
        else warnings_.push_back({lineNo,
            "Unknown grv_version '" + val + "', expected 92|99"});
    }
    // [kbin precision preset]
    else if (key == "precision") {
        std::string vl = toLower(val);
        if      (vl == "custom"      || vl == "0") p.precision = PRECISION_CUSTOM;
        else if (vl == "fast"        || vl == "1") p.precision = PRECISION_FAST;
        else if (vl == "normal"      || vl == "2") p.precision = PRECISION_NORMAL;
        else if (vl == "publication" || vl == "3") p.precision = PRECISION_PUBLICATION;
        else warnings_.push_back({lineNo,
            "Unknown precision '" + val + "', expected custom|fast|normal|publication"});
    }
    // [inclusive mode shortcut]
    else if (key == "inclusive_mode" || key == "inclusive") {
        std::string vl = toLower(val);
        p.inclusive_mode = (vl == "true" || vl == "1" || vl == "yes");
    }
    // [verbosity]
    else if (key == "verbose") {
        std::string vl = toLower(val);
        p.verbose = !(vl == "false" || vl == "0" || vl == "no");
    }
    // [angular cut in degrees -- aliases for cosph_min/max]
    else if (key == "angle_min_deg") {
        double deg = std::stod(val);
        p.cosph_max = std::cos(deg * PI / 180.0);  // min angle -> max cosine
    }
    else if (key == "angle_max_deg") {
        double deg = std::stod(val);
        p.cosph_min = std::cos(deg * PI / 180.0);  // max angle -> min cosine
    }
    // [Run mode]
    else if (key == "mode" || key == "scan_mode") {
        std::string vl = toLower(val);
        if (vl == "x" || vl == "theta" || vl == "kbin" || vl == "all")
            p.scanMode = vl;
        else warnings_.push_back({lineNo,
            "Unknown scan_mode '" + val + "', expected x/theta/kbin/all"});
    }

    // [Integration resolution]
    else if (key == "n_theta_integration" ||
             key == "nth_int")          asInt(p.nth_int);
    else if (key == "n_x_integration" ||
             key == "nx_int")           asInt(p.nx_int);
    // kbin-scan resolution
    else if (key == "ny" || key == "n_y")   asInt(p.ny);
    else if (key == "nkt" || key == "n_kt") asInt(p.nkT);
    // kbin-scan angular acceptance
    else if (key == "cosph_min")            asDouble(p.cosph_min);
    else if (key == "cosph_max")            asDouble(p.cosph_max);

    // unknown key
    else {
        warnings_.push_back({lineNo,
            "Unknown key '" + key + "' -- ignored"});
    }
}

// ============================================================
//  validate  -- physical consistency checks
// ============================================================
std::vector<ValidationError> InputReader::validate(const RunParams &p)
{
    std::vector<ValidationError> errs;

    auto err = [&](const std::string &field, const std::string &msg) {
        errs.push_back({field, msg});
    };

    if (p.E <= 0.0)
        err("E", "Beam energy must be positive");
    if (p.theta_e <= 0.0 || p.theta_e >= 180.0)
        err("theta_e", "Electron scattering angle must be in (0, 180) degrees");
    if (p.L <= 0.0)
        err("Lambda", "Lambda cut-off L must be positive");
    if (p.Ld <= 0.0)
        err("Lambda_Delta", "Lambda cut-off Ld must be positive");
    if (p.xmin < 0.0 || p.xmin >= 1.0)
        err("xmin", "xmin must be in [0, 1)");
    if (p.xmax <= p.xmin || p.xmax > 1.0)
        err("xmax", "xmax must satisfy xmin < xmax <= 1");
    if (p.nx < 2)
        err("nx", "nx must be at least 2");
    if (p.thmin_deg < 0.0 || p.thmin_deg >= 180.0)
        err("theta_min", "theta_min must be in [0, 180) degrees");
    if (p.thmax_deg <= p.thmin_deg || p.thmax_deg > 180.0)
        err("theta_max", "theta_max must satisfy theta_min < theta_max <= 180");
    if (p.xmin_th < 0.0 || p.xmin_th >= 1.0)
        err("x_min_th", "x_min_th must be in [0, 1)");
    if (p.xmax_th <= p.xmin_th || p.xmax_th > 1.0)
        err("x_max_th", "x_max_th must satisfy x_min_th < x_max_th <= 1");
    if (p.kmin_th < 0.0)
        err("k_min_th", "k_min_th must be >= 0");
    if (p.kmax_th <= p.kmin_th)
        err("k_max_th", "k_max_th must be > k_min_th");
    if (p.nth < 2)
        err("nth", "nth must be at least 2");
    if (p.kmin < 0.0)
        err("kmin", "kmin must be >= 0");
    if (p.kmax <= p.kmin)
        err("kmax", "kmax must be > kmin");
    if (p.nk < 2)
        err("nk", "nk must be at least 2");
    if (p.pH < 0.0)
        err("pH", "Hadron momentum pH must be >= 0");
    if (p.alpha1_deg >= p.alpha2_deg)
        err("alpha1/alpha2", "alpha1 must be < alpha2");
    if ((int)p.typ < 0 || (int)p.typ > 6)
        err("form_factor", "form_factor must be in [0,6]");
    if ((int)p.flag < 0 || (int)p.flag > 2)
        err("contribution", "contribution flag must be 0, 1, or 2");
    if ((int)p.dis < 0 || (int)p.dis > 1)
        err("dissociation", "dissociation channel must be 0 or 1");
    if (p.scanMode != "x" && p.scanMode != "theta" &&
        p.scanMode != "kbin" && p.scanMode != "all")
        err("scan_mode", "scan_mode must be 'x', 'theta', 'kbin', or 'all'");
    if (p.cosph_min < 0.0 || p.cosph_min > 1.0)
        err("cosph_min", "cosph_min must be in [0, 1]");
    if (p.cosph_max <= p.cosph_min || p.cosph_max > 1.0)
        err("cosph_max", "cosph_max must satisfy cosph_min < cosph_max <= 1");
    if (p.ny < 2)
        err("ny", "ny must be at least 2");
    if (p.nkT < 2)
        err("nkT", "nkT must be at least 2");

    return errs;
}

// ============================================================
//  printDiagnostics
// ============================================================
void InputReader::printDiagnostics() const
{
    for (auto &w : warnings_)
        std::cerr << "[INPUT WARNING line " << w.line << "] " << w.message << "\n";
    for (auto &e : errors_)
        std::cerr << "[INPUT ERROR field '" << e.field << "'] " << e.message << "\n";
}

// ============================================================
//  writeTemplate  -- write a documented example input file
// ============================================================
void InputReader::writeTemplate(const std::string &filename,
                                const RunParams   &d)
{
    std::ofstream f(filename);
    if (!f) throw std::runtime_error(
        "InputReader::writeTemplate: cannot open '" + filename + "'");

    f <<
"# ============================================================\n"
"#  PionCloud -- Input Parameter File\n"
"#\n"
"#  Syntax:  key = value   (case-insensitive keys)\n"
"#  Lines starting with '#' or ';' are comments.\n"
"#  Section headers [Section] are allowed but ignored.\n"
"#  All parameters are optional; defaults match the original\n"
"#  Fortran programs (3Var_x.f / 3Var_theta.f).\n"
"# ============================================================\n"
"\n"
"[Run Control]\n"
"# Which scan(s) to perform: x | theta | all\n"
"scan_mode        = " << d.scanMode << "\n"
"\n"
"[Output]\n"
"# Output format: none | ascii | root | all\n"
"output_mode      = " << [&]() -> std::string {
    switch (d.outputMode) {
    case OUT_NONE:  return "none";
    case OUT_ASCII: return "ascii";
    case OUT_ROOT:  return "root";
    default:        return "all";
    }
}() << "\n"
"# Base name for all output files (no extension)\n"
"output_name      = " << d.outBaseName << "\n"
"\n"
"[Beam Kinematics]\n"
"# Incident electron beam energy [GeV]\n"
"beam_energy      = " << d.E << "\n"
"# Electron scattering angle [degrees]\n"
"electron_angle   = " << d.theta_e << "\n"
"\n"
"[Hadron Kinematics]  # used by subroutine / f2_pi_sub mode\n"
"# Produced hadron (proton) momentum [GeV]\n"
"hadron_momentum  = " << d.pH << "\n"
"# Hadron production angle lower bound [degrees]\n"
"hadron_angle_min = " << d.alpha1_deg << "\n"
"# Hadron production angle upper bound [degrees]\n"
"hadron_angle_max = " << d.alpha2_deg << "\n"
"\n"
"[Physics - Form Factor]\n"
"# Form factor type (integer code or name):\n"
"#   0 / monopole\n"
"#   1 / dipole\n"
"#   2 / exp_s        (s-dependent exponential, HSS convention)\n"
"#   3 / cov_dipole   (covariant dipole)\n"
"#   4 / dipole_sch   (dipole, s-channel Lambda exchange)\n"
"#   5 / exp_t        (t-dependent exponential, PRD93,054011)\n"
"#   6 / pauli_villar (Pauli-Villar regularisation)\n"
"form_factor      = " << (int)d.typ << "\n"
"\n"
"[Physics - Meson-Baryon Channel]\n"
"# Contribution to include (integer code or name):\n"
"#   0 / pion      -- pi-N     (J = 0 + 1/2)\n"
"#   1 / rho       -- rho-N    (J = 1 + 1/2)\n"
"#   2 / pi_delta  -- pi-Delta (J = 0 + 3/2)\n"
"contribution     = " << (int)d.flag << "\n"
"# Dissociation channel:\n"
"#   0 / charge_exchange  (e.g. N -> pi^- + p)\n"
"#   1 / neutral          (e.g. N -> pi^0 + N)\n"
"dissociation     = " << (int)d.dis << "\n"
"\n"
"[Physics - Cutoff Parameters]\n"
"# Renormalisation cut-off Lambda for pi/rho-N vertex [GeV]\n"
"#   HSS cent.: 1.56  |  HSS+: 1.63  |  HSS-: 1.48\n"
"#   s-dep exp (PRD93,054011): 1.33\n"
"#   t-dep exp (PRD93,054011): 0.85\n"
"#   Pauli-Villar: 0.27\n"
"lambda           = " << d.L << "\n"
"# Renormalisation cut-off Lambda for Delta vertex [GeV]\n"
"#   HSS cent.: 1.39  |  HSS+: 1.46  |  HSS-: 1.32\n"
"lambda_delta     = " << d.Ld << "\n"
"\n"
"[x-Scan Settings]  # active when scan_mode = x or all\n"
"# Minimum Bjorken-x\n"
"x_min            = " << d.xmin << "\n"
"# Maximum Bjorken-x\n"
"x_max            = " << d.xmax << "\n"
"# Number of x output points\n"
"nx               = " << d.nx << "\n"
"# Steps in the theta_h integration (internal, not output points)\n"
"n_theta_integration = " << d.nth_int << "\n"
"\n"
"[theta-Scan Settings]  # active when scan_mode = theta or all\n"
"# Hadron production angle lower bound [degrees]\n"
"theta_min        = " << d.thmin_deg << "\n"
"# Hadron production angle upper bound [degrees]\n"
"theta_max        = " << d.thmax_deg << "\n"
"# Number of theta output points\n"
"nth              = " << d.nth << "\n"
"# Steps in the x integration (internal, not output points)\n"
"n_x_integration  = " << d.nx_int << "\n"
"# x integration bounds for theta-scan [dimensionless]\n"
"# Left panel (Fig.47, full range):       x_min_th=~0,    x_max_th=0.6\n"
"# Right panel (Fig.47, constrained):     x_min_th=0.05,  x_max_th=0.6\n"
"x_min_th         = " << d.xmin_th << "\n"
"x_max_th         = " << d.xmax_th << "\n"
"# |k| integration bounds for theta-scan [GeV]\n"
"# Left panel (Fig.47, full range):       k_min_th=~0,    k_max_th=0.5\n"
"# Right panel (Fig.47, constrained):     k_min_th=0.060, k_max_th=0.250\n"
"k_min_th         = " << d.kmin_th << "\n"
"k_max_th         = " << d.kmax_th << "\n"
"\n"
"[|k| Integration]  # used in both scan modes\n"
"# x-scan:     k integrated over [k_min, k_max] (BoNuS tagging window)\n"
"# theta-scan: k integrated over [k_min_th, k_max_th] (see above)\n"
"# kbin-scan:  k_min/k_max define the |k| bin window\n"
"# Lower bound of |vec{k}| for the x-scan / kbin-scan [GeV]\n"
"#   BoNuS typical bins: 0.060, 0.080, 0.100, 0.130, 0.160 GeV\n"
"k_min            = " << d.kmin << "\n"
"# Upper bound of |vec{k}| for the x-scan / kbin-scan [GeV]\n"
"k_max            = " << d.kmax << "\n"
"# Number of |k| integration steps (used in x-scan and theta-scan)\n"
"nk               = " << d.nk << "\n"
"\n"
"[kbin-Scan Settings]  # active when scan_mode = kbin\n"
"# Mirrors f2pi_sub.f90: integrates (y, kT) at each x\n"
"# applying |k| window [k_min, k_max] and angular acceptance.\n"
"# k_min / k_max above define the |k| bin (e.g. 60-100, 100-200 MeV)\n"
"#\n"
"# Proton polar angle acceptance: cosph_min <= cos(phi_k) <= cosph_max\n"
"#   cos(70°) = 0.342  (upper angle boundary)\n"
"#   cos(30°) = 0.866  (lower angle boundary)\n"
"cosph_min        = " << d.cosph_min << "\n"
"cosph_max        = " << d.cosph_max << "\n"
"# Steps in the y integration (Fortran used 100)\n"
"ny               = " << d.ny << "\n"
"# Steps in the kT integration (Fortran used 50000; use 1000+ for accuracy)\n"
"nkt              = " << d.nkT << "\n";

    std::cout << "[InputReader] Template written to '" << filename << "'\n";
}

} // namespace PionCloud
