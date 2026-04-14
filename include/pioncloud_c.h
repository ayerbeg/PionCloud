/*
 * pioncloud_c.h
 * PionCloud -- C interface for generator / Python / Fortran interoperability
 *
 * This header exposes a flat C API over the C++ library so that programs
 * written in C, Fortran, or Python (via ctypes) can call PionCloud without
 * requiring C++ compilation on the caller side.
 *
 * USAGE FROM C
 * ------------
 *   #include "pioncloud_c.h"
 *   pioncloud_init_default();
 *   double f2pi = pioncloud_eval_f2pi(0.1, 2.0);
 *   double f2n  = pioncloud_eval_f2n (0.1, 2.0);
 *
 * USAGE FROM PYTHON (ctypes)
 * --------------------------
 *   import ctypes, math
 *   lib = ctypes.CDLL("./libPionCloudLib.so")
 *   lib.pioncloud_init_default.restype = None
 *   lib.pioncloud_eval_f2pi.restype    = ctypes.c_double
 *   lib.pioncloud_eval_f2pi.argtypes   = [ctypes.c_double, ctypes.c_double]
 *   lib.pioncloud_init_default()
 *   f2pi = lib.pioncloud_eval_f2pi(0.1, 2.0)
 *
 * USAGE FROM FORTRAN
 * ------------------
 *   interface
 *     subroutine pioncloud_init_default() bind(C)
 *     end subroutine
 *     real(c_double) function pioncloud_eval_f2pi(x, q2) bind(C)
 *       use iso_c_binding
 *       real(c_double), value :: x, q2
 *     end function
 *   end interface
 *
 * THREAD SAFETY WARNING
 * ---------------------
 * CTEQ6 uses Fortran COMMON blocks and SAVE statements.  It is NOT
 * thread-safe.  Do not call these functions from multiple threads
 * simultaneously.  If thread-safe evaluation is needed, serialize
 * calls with a mutex on the calling side.
 *
 * INITIALISATION
 * --------------
 * Call pioncloud_init() or pioncloud_init_default() ONCE before any
 * eval() calls.  Calling init() again with the same parameters is
 * safe (CTEQ6 guards against redundant re-loading).
 */

#ifndef PIONCLOUD_C_H
#define PIONCLOUD_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------
 * pioncloud_params_t
 * Flat parameter struct for the C interface.
 * Mirrors the most commonly needed fields from RunParams.
 * Fields not listed here use the C++ defaults.
 * ---------------------------------------------------------------- */
typedef struct {
    double E;             /* beam energy [GeV]          default: 11.0  */
    double theta_e;       /* electron angle [deg]        default: 35.0  */
    double L;             /* Lambda pi/rho-N [GeV]       default: 1.33  */
    double Ld;            /* Lambda Delta [GeV]          default: 1.39  */
    int    flag;          /* 0=piN, 1=rhoN, 2=piDelta    default: 0     */
    int    typ;           /* form factor 0..6            default: 6     */
    int    dis;           /* 0=charge_exch, 1=neutral    default: 0     */
    int    nucleon;       /* 0=neutron, 1=proton          default: 0     */
    int    grv_version;   /* 92 or 99                    default: 92    */
    double kmin;          /* |k| bin lower [GeV]         default: 0.050 */
    double kmax;          /* |k| bin upper [GeV]         default: 0.100 */
    double cosph_min;     /* cos(phi_k) min              default: 0.342 */
    double cosph_max;     /* cos(phi_k) max              default: 0.866 */
    int    ny;            /* y integration steps         default: 100   */
    int    nkT;           /* kT integration steps        default: 50000 */
    const char *cteq6_table_dir; /* path to .tbl files, NULL=cwd       */
} pioncloud_params_t;

/* ----------------------------------------------------------------
 * pioncloud_init_default
 * Initialise with hardcoded defaults matching the original Fortran:
 *   E=11, theta_e=35, L=1.56, typ=2 (exp_s), flag=0, dis=0
 *   kmin=0.060, kmax=0.250, cosph=[0.342,0.866], ny=100, nkT=5000
 * CTEQ6 table files are read from the current working directory.
 * ---------------------------------------------------------------- */
void pioncloud_init_default(void);

/* ----------------------------------------------------------------
 * pioncloud_init
 * Initialise with user-supplied parameters.
 * ---------------------------------------------------------------- */
void pioncloud_init(const pioncloud_params_t *p);

/* ----------------------------------------------------------------
 * pioncloud_eval_f2pi
 * Returns F2^pi(x, Q2) using the current configuration.
 * Q2 <= 0: computed internally from (E, x, theta_e).
 * ---------------------------------------------------------------- */
double pioncloud_eval_f2pi(double x, double Q2);

/* ----------------------------------------------------------------
 * pioncloud_eval_f2n
 * Returns F2^n(x, Q2) from CTEQ6M via charge symmetry.
 * Returns 0.0 if built without CTEQ6.
 * Q2 <= 0: computed internally.
 * ---------------------------------------------------------------- */
double pioncloud_eval_f2n(double x, double Q2);

/* ----------------------------------------------------------------
 * pioncloud_eval_ratio
 * Returns F2^pi / F2^n at (x, Q2).
 * Returns 0.0 if F2n == 0 or CTEQ6 unavailable.
 * ---------------------------------------------------------------- */
double pioncloud_eval_ratio(double x, double Q2);

/* ----------------------------------------------------------------
 * pioncloud_version
 * Returns a static string describing the build configuration.
 * ---------------------------------------------------------------- */
const char *pioncloud_version(void);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* PIONCLOUD_C_H */
