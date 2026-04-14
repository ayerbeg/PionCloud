// ============================================================
//  pioncloud_c.cc
//  PionCloud C interface implementation
//  Thin wrappers around Calculator and CTEQ6.
// ============================================================

#include "pioncloud_c.h"
#include "PhysicsParams.hh"
#include "Calculator.hh"
#include "Kinematics.hh"

#ifdef HAVE_CTEQ6
#include "CTEQ6.hh"
#endif

#include <cstdlib>
#include <cstring>
#include <string>

// Global Calculator instance used by the C interface.
// Initialised by pioncloud_init() / pioncloud_init_default().
// Not thread-safe (same limitation as CTEQ6).
static PionCloud::Calculator *g_calc = nullptr;
static PionCloud::RunParams   g_params;

// ----------------------------------------------------------------
//  internal: build RunParams from C struct
// ----------------------------------------------------------------
static PionCloud::RunParams makeParams(const pioncloud_params_t *p)
{
    PionCloud::RunParams rp;
    rp.E          = p->E;
    rp.theta_e    = p->theta_e;
    rp.L          = p->L;
    rp.Ld         = p->Ld;
    rp.flag       = static_cast<PionCloud::ContributionFlag>(p->flag);
    rp.typ        = static_cast<PionCloud::FormFactorType>(p->typ);
    rp.dis        = static_cast<PionCloud::DisChannel>(p->dis);
    rp.nucleon    = static_cast<PionCloud::NucleonType>(p->nucleon);
    rp.grv_version= (p->grv_version == 99)
                  ? PionCloud::GRV_99 : PionCloud::GRV_92;
    rp.kmin       = p->kmin;
    rp.kmax       = p->kmax;
    rp.cosph_min  = p->cosph_min;
    rp.cosph_max  = p->cosph_max;
    rp.ny         = p->ny;
    rp.nkT        = p->nkT;
    rp.verbose    = false;   // silent by default in C interface
    rp.outputMode = PionCloud::OUT_NONE;
    if (p->cteq6_table_dir && p->cteq6_table_dir[0] != '\0')
        rp.cteq6_table_dir = std::string(p->cteq6_table_dir);
    return rp;
}

extern "C" {

void pioncloud_init_default(void)
{
    // Defaults matching the original Fortran BoNuS configuration
    PionCloud::RunParams rp;
    rp.E          = 11.0;
    rp.theta_e    = 35.0;
    rp.L          = 1.56;
    rp.Ld         = 1.39;
    rp.flag       = PionCloud::FLAG_PION;
    rp.typ        = PionCloud::FF_EXP_S;
    rp.dis        = PionCloud::DIS_CHARGE_EXCHANGE;
    rp.nucleon    = PionCloud::NUCLEON_NEUTRON;
    rp.grv_version= PionCloud::GRV_92;
    rp.kmin       = 0.060;
    rp.kmax       = 0.250;
    rp.cosph_min  = 0.342;
    rp.cosph_max  = 0.866;
    rp.ny         = 100;
    rp.nkT        = 5000;
    rp.verbose    = false;
    rp.outputMode = PionCloud::OUT_NONE;

    delete g_calc;
    g_params = rp;
    g_calc   = new PionCloud::Calculator(rp);
}

void pioncloud_init(const pioncloud_params_t *p)
{
    if (!p) { pioncloud_init_default(); return; }
    PionCloud::RunParams rp = makeParams(p);
    delete g_calc;
    g_params = rp;
    g_calc   = new PionCloud::Calculator(rp);
}

double pioncloud_eval_f2pi(double x, double Q2)
{
    if (!g_calc) pioncloud_init_default();
    return g_calc->eval(x, Q2);
}

double pioncloud_eval_f2n(double x, double Q2)
{
#ifdef HAVE_CTEQ6
    if (!g_calc) pioncloud_init_default();
    if (Q2 <= 0.0)
        Q2 = PionCloud::calcQ2(g_params.E, x, g_params.theta_e);
    if (Q2 <= 0.0 || x <= 0.0 || x >= 1.0) return 0.0;
    return (g_params.nucleon == PionCloud::NUCLEON_PROTON)
         ? PionCloud::CTEQ6::F2proton (x, Q2, false)
         : PionCloud::CTEQ6::F2neutron(x, Q2, false);
#else
    (void)x; (void)Q2;
    return 0.0;
#endif
}

double pioncloud_eval_ratio(double x, double Q2)
{
    double f2pi = pioncloud_eval_f2pi(x, Q2);
    double f2n  = pioncloud_eval_f2n (x, Q2);
    return (f2n > 0.0) ? f2pi / f2n : 0.0;
}

const char *pioncloud_version(void)
{
    return "PionCloud 1.3"
#ifdef HAVE_CTEQ6
        " +CTEQ6"
#endif
#ifdef HAVE_ROOT
        " +ROOT"
#endif
        ;
}

}  // extern "C"
