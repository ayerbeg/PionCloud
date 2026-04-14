// ============================================================
//  main.cc
//  Pion Cloud Model -- entry point
//
//  PARAMETER PRIORITY (highest wins):
//    1. Command-line flags   (--L, --flag, etc.)
//    2. Input file           (--input myrun.in)
//    3. Built-in defaults    (defined in RunParams{})
//
//  USAGE
//  -----
//    # Normal run from input file
//    ./PionCloud --input myrun.in
//
//    # Generate a template input file and exit
//    ./PionCloud --write-template myrun.in
//
//    # Input file + override one parameter on the command line
//    ./PionCloud --input myrun.in --flag 1 --out MyRhoRun
//
//    # Pure command-line, no input file (all defaults)
//    ./PionCloud --mode x --output ascii --flag 0 --typ 2 --L 1.56
//
//    # List all options
//    ./PionCloud --help
// ============================================================

#include "PhysicsParams.hh"
#include "InputReader.hh"
#include "Calculator.hh"
#include "OutputWriter.hh"

#include <iostream>
#include <string>
#include <cstring>
#include <stdexcept>

using namespace PionCloud;

// ----------------------------------------------------------
//  printUsage
// ----------------------------------------------------------
static void printVersion()
{
    std::cout <<
"PionCloud v1.3 -- Pion Cloud Model Structure Function Calculator\n"
"  C++ port of 3Var_x.f, 3Var_theta.f, f2_pi_sub.f90\n"
"  Build date: " __DATE__ " " __TIME__ "\n"
"  Features compiled in: "
#ifdef HAVE_CTEQ6
"CTEQ6=yes "
#else
"CTEQ6=no "
#endif
#ifdef HAVE_ROOT
"ROOT=yes\n"
#else
"ROOT=no\n"
#endif
;
}

static void printUsage(const char *prog)
{
    std::cout <<
"Usage: " << prog << " [OPTIONS]\n"
"\n"
"  --help / -h          This message\n"
"  --version / -v       Print version and build info\n"
"\n"
"Input file:\n"
"  --input  <file>              Read parameters from INI-style file\n"
"  --write-template <file>      Write a documented template and exit\n"
"\n"
"CLI flags OVERRIDE any value from the input file.\n"
"\n"
"Run control:\n"
"  --mode x|theta|kbin|all      Scan mode (default: all)\n"
"  --verbose 0|1                Print per-point progress (default: 1)\n"
"\n"
"Output:\n"
"  --output ascii|root|all|none Output format (default: all)\n"
"  --out <basename>             Output file base name (default: PionCloud)\n"
"\n"
"Physics:\n"
"  --flag  0|1|2        Contribution: 0=pi-N, 1=rho-N, 2=pi-Delta\n"
"  --dis   0|1          Dissociation: 0=charge exchange, 1=neutral\n"
"  --typ   0..6         Form factor: 0=monopole,1=dipole,2=exp(s),3=cov.dip,\n"
"                                    4=dip(s-ch),5=exp(t),6=Pauli-Villar\n"
"  --L   <GeV>          Lambda for pi/rho-N vertex\n"
"  --Ld  <GeV>          Lambda for Delta vertex\n"
"  --E   <GeV>          Beam energy\n"
"  --theta_e <deg>      Electron scattering angle\n"
"  --nucleon n|p        Nucleon for F2 ratio: neutron (default) or proton\n"
"  --grv 92|99          GRV pion PDF version (default: 92)\n"
"\n"
"x-scan:\n"
"  --xmin <val>         Minimum Bjorken-x\n"
"  --xmax <val>         Maximum Bjorken-x\n"
"  --nx   <n>           Number of output x points\n"
"  --nth-int <n>        Internal theta_h integration steps\n"
"\n"
"theta-scan:\n"
"  --thmin <deg>        theta_h lower bound\n"
"  --thmax <deg>        theta_h upper bound\n"
"  --nth   <n>          Number of output theta points\n"
"  --nx-int  <n>        Internal x integration steps\n"
"  --xmin-th --xmax-th  Internal x bounds for theta-scan\n"
"  --kmin-th --kmax-th  Internal |k| bounds for theta-scan [GeV]\n"
"\n"
"|k| integration / kbin-scan:\n"
"  --kmin <GeV>         Lower |k| bound (x-scan tagging / kbin window)\n"
"  --kmax <GeV>         Upper |k| bound\n"
"  --nk   <n>           |k| steps (x/theta-scan)\n"
"  --ny   <n>           y integration steps (kbin)\n"
"  --nkt  <n>           kT integration steps (kbin; 1000=fast, 50000=exact)\n"
"  --cosph-min <val>    cos(phi_k) lower bound, default 0.342 = cos(70 deg)\n"
"  --cosph-max <val>    cos(phi_k) upper bound, default 0.866 = cos(30 deg)\n"
"  --precision fast|normal|publication  Set ny+nkt preset\n"
"  --inclusive          Remove all kinematic cuts (fully inclusive)\n\n";
}

// ----------------------------------------------------------
//  CliOverrides -- records which flags were explicitly set
//  on the command line so we only overwrite those entries.
// ----------------------------------------------------------
struct CliOverrides {
    bool  has_mode     = false; std::string mode;
    bool  has_output   = false; int    outputMode = OUT_ALL;
    bool  has_out      = false; std::string outName;
    bool  has_flag     = false; int    flag    = 0;
    bool  has_dis      = false; int    dis     = 0;
    bool  has_typ      = false; int    typ     = 6;
    bool  has_L        = false; double L       = 1.33;
    bool  has_Ld       = false; double Ld      = 1.39;
    bool  has_E        = false; double E       = 11.0;
    bool  has_theta_e  = false; double theta_e = 35.0;
    bool  has_pH       = false; double pH      = 0.325;
    bool  has_alpha1   = false; double alpha1  = 30.0;
    bool  has_alpha2   = false; double alpha2  = 70.0;
    bool  has_xmin     = false; double xmin    = 0.02;
    bool  has_xmax     = false; double xmax    = 0.30;
    bool  has_nx       = false; int    nx      = 100;
    bool  has_nth_int  = false; int    nth_int = 100;
    bool  has_thmin    = false; double thmin   = 0.0;
    bool  has_thmax    = false; double thmax   = 100.0;
    bool  has_nth      = false; int    nth     = 100;
    bool  has_nx_int   = false; int    nx_int  = 100;
    bool  has_xmin_th  = false; double xmin_th = 1.0e-4;
    bool  has_xmax_th  = false; double xmax_th = 0.6;
    bool  has_kmin_th  = false; double kmin_th = 1.0e-4;
    bool  has_kmax_th  = false; double kmax_th = 0.5;
    bool  has_kmin     = false; double kmin    = 0.05;
    bool  has_kmax     = false; double kmax    = 0.10;
    bool  has_nk       = false; int    nk      = 100;
    bool  has_ny       = false; int    ny      = 100;
    bool  has_nkT      = false; int    nkT     = 50000;
    bool  has_cosph_min= false; double cosph_min = 0.342;
    bool  has_cosph_max= false; double cosph_max = 0.866;
    bool  has_verbose  = false; bool   verbose = true;
    bool  has_nucleon  = false; NucleonType nucleon = NUCLEON_NEUTRON;
    bool  has_grv      = false; GRVVersion  grv_version = GRV_92;
    bool  has_precision= false; PrecisionPreset precision = PRECISION_CUSTOM;
    bool  has_inclusive= false;
};

static void applyOverrides(RunParams &p, const CliOverrides &o)
{
    if (o.has_mode)    p.scanMode     = o.mode;
    if (o.has_output)  p.outputMode   = o.outputMode;
    if (o.has_out)     p.outBaseName  = o.outName;
    if (o.has_flag)    p.flag         = (ContributionFlag)o.flag;
    if (o.has_dis)     p.dis          = (DisChannel)o.dis;
    if (o.has_typ)     p.typ          = (FormFactorType)o.typ;
    if (o.has_L)       p.L            = o.L;
    if (o.has_Ld)      p.Ld           = o.Ld;
    if (o.has_E)       p.E            = o.E;
    if (o.has_theta_e) p.theta_e      = o.theta_e;
    if (o.has_pH)      p.pH           = o.pH;
    if (o.has_alpha1)  p.alpha1_deg   = o.alpha1;
    if (o.has_alpha2)  p.alpha2_deg   = o.alpha2;
    if (o.has_xmin)    p.xmin         = o.xmin;
    if (o.has_xmax)    p.xmax         = o.xmax;
    if (o.has_nx)      p.nx           = o.nx;
    if (o.has_nth_int) p.nth_int      = o.nth_int;
    if (o.has_thmin)   p.thmin_deg    = o.thmin;
    if (o.has_thmax)   p.thmax_deg    = o.thmax;
    if (o.has_nth)     p.nth          = o.nth;
    if (o.has_nx_int)  p.nx_int       = o.nx_int;
    if (o.has_xmin_th) p.xmin_th      = o.xmin_th;
    if (o.has_xmax_th) p.xmax_th      = o.xmax_th;
    if (o.has_kmin_th) p.kmin_th      = o.kmin_th;
    if (o.has_kmax_th) p.kmax_th      = o.kmax_th;
    if (o.has_kmin)    p.kmin         = o.kmin;
    if (o.has_kmax)    p.kmax         = o.kmax;
    if (o.has_nk)         p.nk           = o.nk;
    if (o.has_ny)         p.ny           = o.ny;
    if (o.has_nkT)        p.nkT          = o.nkT;
    if (o.has_cosph_min)  p.cosph_min    = o.cosph_min;
    if (o.has_cosph_max)  p.cosph_max    = o.cosph_max;
}

// ----------------------------------------------------------
//  main
// ----------------------------------------------------------
int main(int argc, char *argv[])
{
    std::string inputFile;
    std::string templateFile;
    CliOverrides cli;

    for (int i = 1; i < argc; ++i) {
        auto eq  = [&](const char *s){ return std::strcmp(argv[i], s) == 0; };
        auto nxt = [&]() -> std::string {
            if (i+1 >= argc) throw std::runtime_error(
                std::string("Missing argument for ") + argv[i]);
            return argv[++i];
        };
        auto nxtD = [&]() { return std::stod(nxt()); };
        auto nxtI = [&]() { return std::stoi(nxt()); };

        if      (eq("--help")  || eq("-h"))  { printUsage(argv[0]); return 0; }
        else if (eq("--version") || eq("-v")) { printVersion(); return 0; }
        else if (eq("--input") || eq("-i"))  inputFile    = nxt();
        else if (eq("--write-template"))     templateFile = nxt();
        else if (eq("--mode"))    { cli.has_mode    = true; cli.mode      = nxt();  }
        else if (eq("--output")) {
            cli.has_output = true;
            std::string s = nxt();
            if      (s=="ascii") cli.outputMode = OUT_ASCII;
            else if (s=="root")  cli.outputMode = OUT_ROOT;
            else if (s=="all")   cli.outputMode = OUT_ALL;
            else if (s=="none")  cli.outputMode = OUT_NONE;
            else throw std::runtime_error("Unknown output mode: " + s);
        }
        else if (eq("--out"))      { cli.has_out      = true; cli.outName   = nxt();  }
        else if (eq("--flag"))     { cli.has_flag     = true; cli.flag      = nxtI(); }
        else if (eq("--dis"))      { cli.has_dis      = true; cli.dis       = nxtI(); }
        else if (eq("--typ"))      { cli.has_typ      = true; cli.typ       = nxtI(); }
        else if (eq("--L"))        { cli.has_L        = true; cli.L         = nxtD(); }
        else if (eq("--Ld"))       { cli.has_Ld       = true; cli.Ld        = nxtD(); }
        else if (eq("--E"))        { cli.has_E        = true; cli.E         = nxtD(); }
        else if (eq("--theta_e"))  { cli.has_theta_e  = true; cli.theta_e   = nxtD(); }
        else if (eq("--pH"))       { cli.has_pH       = true; cli.pH        = nxtD(); }
        else if (eq("--alpha1"))   { cli.has_alpha1   = true; cli.alpha1    = nxtD(); }
        else if (eq("--alpha2"))   { cli.has_alpha2   = true; cli.alpha2    = nxtD(); }
        else if (eq("--xmin"))     { cli.has_xmin     = true; cli.xmin      = nxtD(); }
        else if (eq("--xmax"))     { cli.has_xmax     = true; cli.xmax      = nxtD(); }
        else if (eq("--nx"))       { cli.has_nx       = true; cli.nx        = nxtI(); }
        else if (eq("--nth-int"))  { cli.has_nth_int  = true; cli.nth_int   = nxtI(); }
        else if (eq("--thmin"))    { cli.has_thmin    = true; cli.thmin     = nxtD(); }
        else if (eq("--thmax"))    { cli.has_thmax    = true; cli.thmax     = nxtD(); }
        else if (eq("--nth"))      { cli.has_nth      = true; cli.nth       = nxtI(); }
        else if (eq("--nx-int"))   { cli.has_nx_int   = true; cli.nx_int    = nxtI(); }
        else if (eq("--kmin"))     { cli.has_kmin     = true; cli.kmin      = nxtD(); }
        else if (eq("--kmax"))     { cli.has_kmax     = true; cli.kmax      = nxtD(); }
        else if (eq("--nk"))       { cli.has_nk       = true; cli.nk        = nxtI(); }
        else if (eq("--ny"))       { cli.has_ny       = true; cli.ny        = nxtI(); }
        else if (eq("--nkt"))      { cli.has_nkT      = true; cli.nkT       = nxtI(); }
        else if (eq("--cosph-min")){ cli.has_cosph_min= true; cli.cosph_min = nxtD(); }
        else if (eq("--cosph-max")){ cli.has_cosph_max= true; cli.cosph_max = nxtD(); }
        else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            printUsage(argv[0]); return 1;
        }
    }

    // ---- --write-template mode ----
    if (!templateFile.empty()) {
        RunParams defaults;
        applyOverrides(defaults, cli);
        InputReader::writeTemplate(templateFile, defaults);
        return 0;
    }

    // ---- build RunParams: defaults -> file -> CLI ----
    RunParams params;  // starts with all defaults

    if (!inputFile.empty()) {
        InputReader reader;
        try {
            reader.readInto(inputFile, params);
        } catch (const std::exception &e) {
            std::cerr << "Fatal: " << e.what() << "\n";
            return 1;
        }
        reader.printDiagnostics();
        if (reader.hasErrors()) {
            std::cerr << "Input file has validation errors -- aborting.\n";
            return 1;
        }
    }

    applyOverrides(params, cli);  // CLI flags win over file

    // Final validation (catches bad CLI values too)
    auto errs = InputReader::validate(params);
    if (!errs.empty()) {
        for (auto &e : errs)
            std::cerr << "[VALIDATION ERROR] " << e.field
                      << ": " << e.message << "\n";
        return 1;
    }

    // ---- print active configuration ----
    std::string cfgSrc = inputFile.empty()
        ? "(defaults + CLI)"
        : "(file: " + inputFile + " + CLI overrides)";

    // ---- helpers to print enum names ----
    auto ffNameStr = [](FormFactorType t) -> std::string {
        switch (t) {
        case FF_MONOPOLE:     return "monopole (0)";
        case FF_DIPOLE:       return "dipole (1)";
        case FF_EXP_S:        return "exp_s (2)";
        case FF_COV_DIPOLE:   return "cov_dipole (3)";
        case FF_DIPOLE_S_CH:  return "dipole_sch (4)";
        case FF_EXP_T:        return "exp_t (5)";
        case FF_PAULI_VILLAR: return "pauli_villar (6)";
        default:              return "unknown";
        }
    };
    auto flagNameStr = [](ContributionFlag f) -> std::string {
        switch (f) {
        case FLAG_PION:      return "pi-N (0)";
        case FLAG_RHO:       return "rho-N (1)";
        case FLAG_PI_DELTA:  return "pi-Delta (2)";
        default:             return "unknown";
        }
    };
    auto disNameStr = [](DisChannel d) -> std::string {
        switch (d) {
        case DIS_CHARGE_EXCHANGE: return "charge exchange (0)";
        case DIS_NEUTRAL:         return "neutral (1)";
        default:                  return "unknown";
        }
    };
    auto outNameStr = [](int m) -> std::string {
        switch (m) {
        case OUT_NONE:  return "none";
        case OUT_ASCII: return "ascii";
        case OUT_ROOT:  return "root";
        default:        return "all";
        }
    };

    std::cout <<
"\n╔══════════════════════════════════════════════╗\n"
"║  Pion Cloud Model -- Structure Function Calc ║\n"
"╚══════════════════════════════════════════════╝\n"
"\nActive configuration " << cfgSrc << ":\n"
"  scan_mode      : " << params.scanMode    << "\n"
"  E              : " << params.E           << " GeV\n"
"  theta_e        : " << params.theta_e     << " deg\n"
"  pH             : " << params.pH          << " GeV\n"
"  alpha1/2       : " << params.alpha1_deg  << " / "
                      << params.alpha2_deg  << " deg\n"
"  Lambda (L)     : " << params.L           << " GeV\n"
"  Lambda (Ld)    : " << params.Ld          << " GeV\n"
"  Form factor    : " << ffNameStr(params.typ)    << "\n"
"  Contribution   : " << flagNameStr(params.flag) << "\n"
"  Dis. channel   : " << disNameStr(params.dis)   << "\n"
"  x range        : [" << params.xmin << ", " << params.xmax
                       << "]  nx=" << params.nx
                       << "  nth_int=" << params.nth_int << "\n"
"  theta range    : [" << params.thmin_deg << ", " << params.thmax_deg
                       << "] deg  nth=" << params.nth
                       << "  nx_int=" << params.nx_int << "\n"
"  theta x range  : [" << params.xmin_th << ", " << params.xmax_th << "]\n"
"  theta k range  : [" << params.kmin_th << ", " << params.kmax_th << "] GeV\n"
"  |k| range      : [" << params.kmin << ", " << params.kmax
                       << "] GeV  nk=" << params.nk << "\n"
"  kbin settings  : ny=" << params.ny
                       << "  nkT=" << params.nkT
                       << "  cosph=[" << params.cosph_min
                       << ", " << params.cosph_max << "]\n"
"  output_mode    : " << outNameStr(params.outputMode) << "\n"
"  output_name    : " << params.outBaseName << "\n"
"  F2n / ratio    : "
#ifdef HAVE_CTEQ6
  "enabled (CTEQ6M, iset=1)\n"
#else
  "disabled (rebuild with -DUSE_CTEQ6=ON)\n"
#endif
"\n";

    // ---- run calculations ----
    Calculator calc(params);
    ScanResults results;

    try {
        if (params.scanMode == "x" || params.scanMode == "all") {
            std::cout << "=== Running x-scan ===\n";
            calc.runXScan(results);
        }
        if (params.scanMode == "theta" || params.scanMode == "all") {
            std::cout << "\n=== Running theta-scan ===\n";
            calc.runThetaScan(results);
        }
        if (params.scanMode == "kbin" || params.scanMode == "all") {
            std::cout << "\n=== Running kbin-scan ===\n";
            calc.runKBinScan(results);
        }
    } catch (const std::exception &e) {
        std::cerr << "Calculation error: " << e.what() << "\n";
        return 1;
    }

    // ---- write output ----
    OutputWriter writer(params, params.outBaseName);
    writer.write(results);

    std::cout << "\nDone.\n";
    return 0;
}
