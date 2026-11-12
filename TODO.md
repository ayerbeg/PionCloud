# PionCloud -- TODO / Roadmap
# ============================================================
# Items raised during development and code review sessions.
# Priority: HIGH / MEDIUM / LOW
# Status:   open / in-progress / done
# ============================================================

## FINDINGS OF THE 2026-11-12 REVIEW  (details: analyses/code_review/REVIEW.md on the laptop)
# --------------------------------------------------------

[HIGH] [done v1.4]  GRV99 DOUBLE-COUNTS THE PION VALENCE
  GRV99 (GRS99) xV is the total valence u_v+dbar_v (int V dx = 1.99),
  GRV92 xV is one flavour (int V dx = 0.99); both enter 5/9(xV+2xS), so
  grv_version=99 gives F2piN x2 (+107% on Fig.37).  Use xV/2 for GRV99 and
  check the sea normalisation vs EPJ C10 (1999) 313.  In f2_pi_sub.f90 the
  GRV99 call was commented out (never validated).  Add a sum-rule test.
  DONE: GRS99 Appendix A implemented (also: strange sea was used as light
  sea; sign of 0.009 s^2).  Sum-rule tests added.  See CHANGELOG.md.

[MEDIUM] [done v1.4]  DEFAULT FORM FACTOR: PAULI-VILLARS WITH L=1.33
  Confirmed with PRD 93 054011 Table I: 1.33 = s-exp fit, PV fit = 0.27.
  Default now exp_s / 1.56 (HSS); old default gave F2 x1.6-1.9.  PV
  implementation = their Eq. (9) (correct).  Range warning added.

[MEDIUM] [done v1.4]  DELTA CHANNELS RETURN 0 FOR dissociation=neutral (silently)
  Now rejected by InputReader::validate and Calculator.
[LOW] [open]  EXPLICIT DELTA CHARGE STATES (p/n target x final state with
  Clebsch-Gordan factors) instead of the Fortran's summed factor 2/3.

[LOW] [done v1.4]  rho-Delta (f_RhoDel) IMPLEMENTED BUT NOT SELECTABLE  -> flag 3
[LOW] [done v1.4]  ODD SIMPSON STEP COUNTS -> Simpson + 3/8 tail, exact for any n
[HIGH] [done v1.4]  CTEST FALSE PASS: Fortran STOP (exit 0) on missing table
  after ~80 tests; now WORKING_DIRECTORY + PASS_REGULAR_EXPRESSION, and the
  Calculator throws before calling SetCtq6.  197 tests really run.
[MEDIUM] [done v1.4]  9 CLI OPTIONS IN --help WERE NEVER PARSED (--grv,
  --nucleon, --precision, --inclusive, --verbose, --x/kmin/max-th)
[LOW] [open]  KBIN: ANALYTIC kT LIMITS + GAUSS QUADRATURE INSTEAD OF
              kT in [0,10] GeV WITH continue-CUTS (accuracy, speed)
[LOW] [open]  ARCHITECTURE: SplittingFunction / Regulator / MesonSF /
              Region interfaces + registry; LHAPDF backend; one kernel
              for all scan modes (see review Part 5)

## PHYSICS OUTPUT
# --------------------------------------------------------

[HIGH] [done]  NUCLEON SELECTOR
  nucleon = neutron|proton in RunParams and input file.
  F2proton() and F2neutron() wired; nucleon selector affects
  CTEQ6 call in all three scan modes and ROOT/ASCII output labels.

[HIGH] [done]  GRV99 PION PDF (UPDATED PARAMETRISATION)
  GRV99() added to GRV.hh, translated from f2_pi_sub.f90 lines 766-823.
  Selected via grv_version = 92|99 in input file.
  Used in integrand() (x/theta-scan) and runKBinScan() / eval().

[MEDIUM] [done]  RATIO F2pi/F2p (PROTON DENOMINATOR)
  Follows from nucleon selector above -- set nucleon=proton.

[MEDIUM] [done]  INCLUSIVE CURVE (NO CUTS) AS A DEDICATED MODE
  inclusive_mode = true in input file removes all kinematic cuts
  (cosph_min=0, cosph_max=1, kmin=0.001, kmax=10.0).

[MEDIUM] [open]  SUM OF CHANNELS (piN + rhoN + piDelta)
  The total pion-cloud contribution F2^(piN)(x) is the sum of the
  three individual channels.  Currently user must run three separate
  scans and add results manually.  A combined scan_mode = all_channels
  that runs all three and stores the sum would be useful.
  Would require new ScanResults field for the sum.

[LOW] [done]  ANGULAR CUT AS DEGREES IN INPUT FILE
  angle_min_deg / angle_max_deg aliases added to InputReader.
  Converts to cosph internally: cosph = cos(deg * pi/180).


## OUTPUT / ROOT
# --------------------------------------------------------

[HIGH] [done]  THETA-SCAN OBJECTS USE FIXED NAMES
  g_thetascan etc. now use the <channel>_<mode> prefix:
  e.g. g_piN_neutral_theta, g_rhoN_neutral_theta.

[MEDIUM] [done]  LOG SCALE NOT SET ON THETA-SCAN CANVAS
  c_<pfx>_theta now calls SetLogy(1).

[MEDIUM] [done]  ASCII XSCAN.DAT COLUMN HEADER IS ALWAYS "F2piK"
  Column header now encodes channel and mode:
  F2piK_piN_kbin, F2piK_piN_neutral_theta, etc.

[HIGH] [done]  X-SCAN/KBIN ASCII FILE NAME AND COLUMN SEPARATION (2026-11)
  The x-scan/kbin table was written as "<base>_pi-N (flag=0)_kbin.dat"
  (spaces/parentheses, label from flagName instead of channelTag), and
  the F2n column (setw 14 for a 14-char %.8e number) was glued to the
  running moment, so the file could not be read with whitespace parsers.
  Now <base>_xscan.dat as documented, labels use channelTag, all fields
  are separated, and the F2N comment no longer claims "charge symmetry"
  for nucleon=proton (params file said CTEQ6L; the code uses CTEQ6M).

[LOW] [open]  RUNNING MOMENT IN THETA-SCAN IS TRAPEZOIDAL
  The theta-scan running moment is trapezoidal while the first
  moment uses Simpson.  The x-scan running moment was updated to
  use Simpson (consistent with F2pi0_x).  Theta-scan still uses
  trapezoidal -- fix to use Simpson for full consistency.


## INTEGRATION / NUMERICS
# --------------------------------------------------------

[HIGH] [done]  KBIN nkT DEFAULT IS 50000 -- VERY SLOW
  Precision preset added: precision = fast|normal|publication
  sets ny and nkT simultaneously:
    fast:        ny=20,  nkT=500
    normal:      ny=50,  nkT=2000
    publication: ny=100, nkT=50000 (matches Fortran, default)
  Documented in --help and template file.

[MEDIUM] [done]  KBIN SCAN x RANGE: SMALL-x Q2 WARNING
  If Q2 <= 0 at any x point (can happen at theta_e=12 deg, x<0.01),
  a warning is printed to stderr and the point is skipped gracefully.


## USABILITY
# --------------------------------------------------------

[MEDIUM] [done]  NO --version FLAG
  --version / -v prints version, build date, and compiled features.

[MEDIUM] [done]  MISSING CLI FLAGS FOR KBIN IN HELP TEXT
  --ny, --nkt, --cosph-min, --cosph-max, --precision, --inclusive,
  --verbose, --nucleon, --grv all listed in printUsage().
  NOTE (v1.4): v1.3 only LISTED them; --precision, --inclusive, --verbose,
  --nucleon, --grv were not parsed.  Implemented in v1.4.

[LOW] [open]  TEMPLATE INPUT FILE DOES NOT INCLUDE KBIN EXAMPLE
  --write-template generates the [kbin-Scan Settings] section with
  all keys documented, but no worked example showing the four k-bin
  curves as a comment block at the top.


## TESTING
# --------------------------------------------------------

[MEDIUM] [done]  NO NUMERICAL REGRESSION TESTS
  Three regression tests pinned against real Fortran output:
    1. runKBinScan: k=[200,300] MeV, 12 points (bin t=2)
    2. runKBinScan: k=[300,400] MeV, 15 points (bin t=3)
    3. runKBinScan: k=[60,250] MeV,   7 points (f2pi_sub, 30-70 deg cut)
  All three from f2_pi_ranges.c / f2_pi.c run on user's machine.
  Still missing:
    - runXScan regression (3Var_x.f): need 3Var_x.f output with
      k=[60,250] MeV -- this is the actual black piN curve in f2piN.png
    - runThetaScan regression (3Var_theta.f): need Fig.47 Fortran output
  Partial (2026-11): docs/validation/ compares x, theta and kbin scans
  with curves extracted from the vector figures of PR12-15-006 (Figs. 5,
  6, 37, 46, 47): agreement ~1% (bulk <= 2-4%).  Not yet run by ctest.

[MEDIUM] [open]  CTEST REGRESSION AGAINST PUBLISHED CURVES
  Turn docs/validation/ref_data into a ctest (tolerance ~3% in the bulk,
  where the curve is within x100 of its maximum).  Needs no ROOT/Python
  if done in C++ with the reference .dat files.

[LOW] [open]  F2n / F2p NORMALISATION VS PUBLISHED ORANGE CURVES
  CTEQ6M F2 (used for the ratio) is 3-7% above the F2n/F2p curves of
  PR12-15-006 Figs. 5, 6, 37 (x-dependent, max near x~0.15).  Not due
  to theta_e.  Find which set/order/Q2 the published curves used.

[LOW] [done]  TEST COVERAGE FOR SCAN_MODE=ALL WITH KBIN
  Documented: scan_mode=all runs kbin last, which overwrites xScan.
  Tests now use explicit scan modes instead of relying on "all"
  for numerical comparisons.  runAll() structural test still present.

[LOW] [open]  CTEQ6 TESTS NOT TESTED WITH CTEQ6 ACTIVE
  Template round-trip and new keys are tested in no-CTEQ6 build only.


## GENERATOR / LIBRARY INTERFACE
# --------------------------------------------------------

[HIGH] [done]  SINGLE-POINT EVALUATION METHOD
  double Calculator::eval(double x, double Q2 = -1.0) const
  Runs kbin (y,kT) integral for one x point. Q2<=0 triggers
  internal calcQ2(). Always silent. Tested against runKBinScan.

[HIGH] [done]  VERBOSITY CONTROL
  RunParams::verbose = true|false (default: true).
  All std::cout in Calculator.cc and OutputWriter.cc guarded.
  Input key: verbose = true|false|0|1|yes|no.

[HIGH] [done]  CTEQ6 TABLE PATH CONFIGURABLE
  RunParams::cteq6_table_dir = "/path/to/tables".
  Calculator constructor chdir()s there before CTEQ6::init()
  and restores the working directory afterwards.

[MEDIUM] [done]  C INTERFACE (extern "C" WRAPPER)
  include/pioncloud_c.h  -- flat C API with usage examples for C,
                            Python (ctypes), and Fortran.
  src/pioncloud_c.cc     -- implementation with global Calculator.
  Functions: pioncloud_init_default(), pioncloud_init(),
             pioncloud_eval_f2pi(), pioncloud_eval_f2n(),
             pioncloud_eval_ratio(), pioncloud_version().

[MEDIUM] [open]  THREAD SAFETY DOCUMENTATION + WARNING
  CTEQ6 uses Fortran COMMON blocks -- not thread-safe.
  pioncloud_c.h has a warning comment. Full mutex wrapper or
  per-thread initialisation not yet implemented.

[LOW] [open]  PYTHON BINDINGS (pybind11 or ctypes)
  C interface is complete; Python ctypes usage example is shown
  in pioncloud_c.h.  A proper pybind11 or ctypes wrapper module
  with pip-installable packaging is not yet implemented.
  Requires pybind11 to be available at build time.


## DOCUMENTATION
# --------------------------------------------------------

[LOW] [done]  README DOES NOT EXPLAIN THE PHYSICS BEHIND EACH MODE
  LaTeX manual (docs/PionCloud_manual.tex + .pdf) covers physics
  background, scan modes with equations, and worked examples.

[LOW] [open]  CTEQ6/README SHOULD MENTION THE PROTON VS NEUTRON OPTION
  Update CTEQ6/README to note that nucleon=proton uses F2proton()
  directly (no charge symmetry), while nucleon=neutron uses
  charge symmetry: u_n=d_p, d_n=u_p.

[LOW] [done]  MANUAL: FIGURES AND CODE FONT (2026-11)
  Manual referred to "Fig. 47", "black curve", "four k-bin curves"
  without showing or citing them.  New section "Worked examples: the
  published TDIS figures" cites PR12-15-006 / PR12-14-010 / Hobbs
  arXiv:1411.2216, gives inputs and overlays PionCloud on the published
  curves (docs/figures, docs/validation).  Corrected: HSS = Holtmann-
  Szczurek-Speth (not "Hobbs, Shanahan, Strikman"); theta_e = 35 deg also
  for kbin (12 deg does not reproduce the figures); inclusive_mode gives
  one charge state (x3 neutral, x3/2 charge exchange for the total);
  pi-Delta uses lambda_delta; nkt=50000 takes ~1 s, not hours.  Listings
  now use Inconsolata, upright comments, columns=fullflexible (the old
  CM-typewriter italic was letter-spaced by columns=fixed).

[LOW] [open]  KSCAN MODE (F2piN vs |k|)
  Needed to reproduce Fig. 45 of PR12-15-006 (|k| dependence for
  Delta_x=[0,0.6], all angles).  Reference curves can be extracted
  from page 66 with docs/validation/extract_curves.py.

[LOW] [open]  TSCAN MODE (F2piN vs -t / Mandelstam)
  Plot F2^(piN) vs -t by binning the (y,kT) integrand over
  t(kT,y) = -kT^2 - mN^2*y^2/(1-y) windows.
  Discussed but not yet designed in detail.
  Open question: fixed x or x-integrated?

# ============================================================
# END OF TODO LIST
# Last updated: November 2026
# Version: v1.4 (in progress)
# Tests: 197 passed, 0 failed (v1.4; all now actually executed -- see CHANGELOG)
# ============================================================
