#pragma once
// ============================================================
//  CTEQ6.hh
//  C++ interface to the CTEQ6 Fortran PDF library
//  (Pumplin et al., JHEP 0207:012, 2002 -- 2007 update v6.52)
//
//  Source file: CTEQ6/Cteq6Pdf-2007.f
//
//  PUBLIC FORTRAN API (the only two symbols used externally):
//    SetCtq6 (Iset)           -- subroutine, initialises PDF tables
//    Ctq6Pdf (Iparton, X, Q)  -- function,   returns xf(x,Q)
//
//  PRECISION
//  ---------
//  Cteq6Pdf-2007.f uses "Implicit REAL (A-H,O-Z)", meaning all
//  real variables are 4-byte single precision.  Correspondingly:
//    - Ctq6Pdf returns a C float (not double)
//    - X and Q arguments must be passed as float* (not double*)
//  This header handles all casts internally; callers always use
//  double on the C++ side.
//
//  GFORTRAN NAME MANGLING
//  ----------------------
//  gfortran appends a trailing underscore to Fortran symbol names:
//    SetCtq6  ->  setctq6_
//    Ctq6Pdf  ->  ctq6pdf_
//  If you use ifort or flang, verify their mangling convention and
//  adjust the extern "C" declarations below accordingly.
//
//  ISET TABLE (sets available with your four .tbl files)
//  -------------------------------------------------------
//  Iset  Name      Scheme  alphaS(Mz)  Table file
//  ----  --------  ------  ----------  -----------
//   1    CTEQ6M    MSbar   0.118       cteq6m.tbl   (NLO central -- used by Fortran)
//   2    CTEQ6D    DIS     0.118       cteq6d.tbl   (NLO DIS scheme)
//   3    CTEQ6L    LO      0.118       cteq6l.tbl   (LO, same alphaS as NLO)
//   4    CTEQ6L1   LO      0.130       cteq6l1.tbl  (LO, alphaS from LO fit)
//
//  The original Fortran programs call SETCTQ6(1) with the comment
//  "CTEQ MS-bar scheme", so Iset=1 (CTEQ6M) is the default here.
//
//  PARTON CODES (Iparton argument)
//  --------------------------------
//   0   gluon
//   1   u        -1  u-bar
//   2   d        -2  d-bar
//   3   s        -3  s-bar
//   4   c        -4  c-bar
//   5   b        -5  b-bar
//
//  TABLE FILES
//  -----------
//  SetCtq6 opens the table file by name in the current working
//  directory (no path prefix).  CMake copies the *.tbl files from
//  CTEQ6/ to the build directory at configure time, so running
//  the executable from the build directory works automatically.
//  If you run from elsewhere, either set the working directory or
//  symlink the .tbl files next to the executable.
//
//  USAGE
//  -----
//    #include "CTEQ6.hh"
//    PionCloud::CTEQ6::init();                   // call once
//    double f2n = PionCloud::CTEQ6::F2neutron(x, Q2);
//    double rat = F2piK / f2n;
// ============================================================

#include <cmath>

// ----------------------------------------------------------
//  Raw Fortran bindings -- not for direct use outside this header
// ----------------------------------------------------------
extern "C" {
    // subroutine SetCtq6(Iset)
    void setctq6_(int *iset);

    // function Ctq6Pdf(Iparton, X, Q)  -- returns REAL (float)
    float ctq6pdf_(int *iparton, float *x, float *q);
}

namespace PionCloud {
namespace CTEQ6 {

// ----------------------------------------------------------
//  PDF set identifiers -- matches Iset values in SetCtq6
// ----------------------------------------------------------
enum Set {
    CTEQ6M  = 1,   // NLO MSbar, alphaS(Mz)=0.118  -- Fortran default
    CTEQ6D  = 2,   // NLO DIS scheme
    CTEQ6L  = 3,   // LO, alphaS(Mz)=0.118
    CTEQ6L1 = 4    // LO, alphaS(Mz)=0.130 (from LO fit)
};

// ----------------------------------------------------------
//  init()
//  Load the PDF grid tables.  Must be called once before any
//  PDF query.  Safe to call multiple times (Fortran guards with
//  a saved Isetold flag and skips re-loading if Iset unchanged).
//
//  Default: CTEQ6M (iset=1), matching the original Fortran:
//    CALL SETCTQ6(1)  ! CTEQ 'MS-bar' SCHEME
// ----------------------------------------------------------
inline void init(Set set = CTEQ6M)
{
    int iset = static_cast<int>(set);
    setctq6_(&iset);
}

// ----------------------------------------------------------
//  xfx()
//  Returns x * f(x, Q) for the requested parton.
//  x        -- Bjorken-x          (valid range: ~1e-6 to 1)
//  Q        -- factorisation scale [GeV]  (valid range: >1.3 GeV)
//  iparton  -- parton code (see table above)
//
//  Casts to/from float internally; the Fortran library is
//  single-precision throughout.
// ----------------------------------------------------------
inline double xfx(int iparton, double x, double Q)
{
    float fx = static_cast<float>(x);
    float fQ = static_cast<float>(Q);
    int   ip = iparton;
    return static_cast<double>(ctq6pdf_(&ip, &fx, &fQ));
}

// ----------------------------------------------------------
//  F2proton()
//  Proton structure function F₂^p(x, Q²) from CTEQ6 PDFs.
//
//  F₂^p = x * [(4/9)(u+ū) + (1/9)(d+d̄) + (1/9)(s+s̄)]
//
//  Strange quarks are optionally included (default: false) to
//  match the original Fortran which omits them.
// ----------------------------------------------------------
inline double F2proton(double x, double Q2, bool include_strange = false)
{
    if (x  <= 0.0 || x  >= 1.0) return 0.0;
    if (Q2 <= 0.0)               return 0.0;
    const double Q = std::sqrt(Q2);

    double u    = xfx( 1, x, Q);
    double ubar = xfx(-1, x, Q);
    double d    = xfx( 2, x, Q);
    double dbar = xfx(-2, x, Q);

    double f2p = x * ( (4.0/9.0)*(u + ubar)
                     + (1.0/9.0)*(d + dbar) );
    if (include_strange) {
        double s    = xfx( 3, x, Q);
        double sbar = xfx(-3, x, Q);
        f2p += x * (1.0/9.0)*(s + sbar);
    }
    return f2p;
}

// ----------------------------------------------------------
//  F2neutron()
//  Neutron structure function F₂ⁿ(x, Q²) via charge symmetry.
//
//  Charge symmetry: u_n = d_p, d_n = u_p, s_n = s_p
//  Therefore:
//    F₂ⁿ = x * [(4/9)(d+d̄) + (1/9)(u+ū)]
//
//  This matches the original Fortran exactly:
//    F2neu = x * ((4/9)*(d_pro+dbar_pro) + (1/9)*(u_pro+ubar_pro))
//  (f2_pi_sub.f90 line 269, 3Var_x.f line 186)
//
//  Note: the factor of 2 in 3Var_x.f (F2neu = 2.*x*...) is a
//  bug corrected in f2_pi_sub.f90 (comment: "corrected from
//  original").  We follow the corrected version: factor = x only.
// ----------------------------------------------------------
inline double F2neutron(double x, double Q2, bool include_strange = false)
{
    if (x  <= 0.0 || x  >= 1.0) return 0.0;
    if (Q2 <= 0.0)               return 0.0;
    const double Q = std::sqrt(Q2);

    double u    = xfx( 1, x, Q);
    double ubar = xfx(-1, x, Q);
    double d    = xfx( 2, x, Q);
    double dbar = xfx(-2, x, Q);

    double f2n = x * ( (4.0/9.0)*(d + dbar)
                     + (1.0/9.0)*(u + ubar) );
    if (include_strange) {
        double s    = xfx( 3, x, Q);
        double sbar = xfx(-3, x, Q);
        f2n += x * (1.0/9.0)*(s + sbar);
    }
    return f2n;
}

} // namespace CTEQ6
} // namespace PionCloud
