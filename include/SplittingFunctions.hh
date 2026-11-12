#pragma once
// ============================================================
//  SplittingFunctions.hh
//  Pion Cloud Model -- meson-baryon splitting functions
//
//  Direct C++ translations of the Fortran functions of T. Hobbs
//  (3Var_x.f, 3Var_theta.f, f2_pi_sub.f90):
//    fypiN    -- pi-N      (N -> pi N),     cut-off L
//    f_rhoN   -- rho-N     (N -> rho N),    cut-off L
//    fypiD    -- pi-Delta  (N -> pi Delta), cut-off Ld
//    f_RhoDel -- rho-Delta (N -> rho Delta),cut-off Ld
//
//  WHAT THEY ARE
//  f_MB(y, kT) is the light-cone (IMF / time-ordered PT) probability
//  density for the nucleon to be found as meson M + baryon B, with the
//  meson carrying light-cone momentum fraction y and transverse momentum
//  kT.  Normalisation:
//      <n>_MB = int_0^1 dy int_0^inf d(kT^2) f_MB(y, kT)
//  i.e. the kT measure is d(kT^2) = 2 kT dkT (this is why the kbin scan
//  integrates 2*kT*f).  The Sullivan contribution to F2 is
//      F2^(MB)(x) = int_x^1 dy int d(kT^2) f_MB(y,kT) F2^M(x/y, Q2).
//  References: Holtmann, Szczurek, Speth NPA 596 (1996) 631;
//  Melnitchouk, Speth, Thomas PRD 59 (1999) 014033; Hobbs, Londergan,
//  Melnitchouk PRD 89 (2014) 074008; McKenney et al. PRD 93 (2016) 054011.
//
//  VERIFIED (v1.4, analyses/code_review/multiplicity): summed over charge
//  states, exp_s Lambda=1.33 gives <n>_piN = 0.224 and Pauli-Villars
//  Lambda=0.27 gives 0.310, vs 0.23 and 0.31 in McKenney et al. Table I.
//  OPEN: <n>_piDelta = 0.027 vs their 0.06 (factor ~2.2) -- see TODO.md.
//
//  VARIABLES
//    S    : invariant mass^2 of the MB state,
//           S = (kT^2 + m_M^2)/y + (kT^2 + m_B^2)/(1-y)
//    t    : meson virtuality, t = -(kT^2 + y^2 mN^2)/(1-y)   (pi-N)
//  FORM FACTORS (formFactor below; L = Lambda):
//    exp_s : exp[(mN^2 - S)/L^2]   (HSS / PR12-15-006 Eq. 14 convention;
//            HSS 1996 use exp[(M^2-S)/2L^2], i.e. their L is smaller by sqrt2)
//    exp_t : exp[(t - m_M^2)/L^2]
//    pauli_villar : sqrt(1 - (t-m_M^2)^2/(t-L^2)^2)  = McKenney Eq. (9)
//    others: see formFactor().
//  COUPLINGS (values as in the Fortran):
//    g_piNN^2/4pi = 14.40 ; g_rhoNN^2/4pi = 0.55 with kappa = f/g = 6.1 ;
//    pi-N-Delta: 0.2237 (x 4pi, with an explicit 1/m_pi^2) ;
//    rho-N-Delta: 20.448 (x 4pi, with an explicit 1/m_rho^2).
//  ISOSPIN (argument dis):
//    pi-N, rho-N : charge_exchange (n -> M^- p, p -> M^+ n) factor 2,
//                  neutral (N -> M^0 N) factor 1; sum over states = 3.
//    pi-Delta, rho-Delta : charge_exchange = Fortran sum over charge
//                  states (factor 2/3); neutral -> coupling 0 (rejected by
//                  Calculator/InputReader since v1.4).
//
//  HOW TO ADD A CHANNEL
//    1. write   inline double f_XY(double y, double kT, double L,
//                                  FormFactorType typ, DisChannel dis)
//       here, returning f per d(kT^2) as above;
//    2. add a ContributionFlag value in PhysicsParams.hh;
//    3. add one case to Calculator::splitting() (Calculator.cc), plus the
//       names in InputReader.cc / OutputWriter.cc / main.cc;
//    4. check <n> against a published value and add a unit test.
//  (Library users can instead inject a model at run time, see Models.hh.)
//
//  All functions are stateless and thread-safe.
//  Arguments:
//    y   -- light-cone momentum fraction carried by the meson
//    kT  -- meson transverse momentum [GeV]
//    L   -- form-factor cut-off Lambda [GeV]
//    typ -- FormFactorType enum
//    dis -- DisChannel enum
// ============================================================

#include "PhysicsParams.hh"
#include <cmath>
#include <stdexcept>

namespace PionCloud {

// ----------------------------------------------------------
//  Form-factor helper
//    Returns FF value for a given invariant mass S, type, L
//    and the "t" variable (needed for some types).
// ----------------------------------------------------------
inline double formFactor(double S, double t,
                         double mMeson, double mBaryon,
                         double L, FormFactorType typ,
                         double sM = 0.0)
{
    switch (typ) {
    case FF_MONOPOLE:
        return (L*L + mN*mN) / (L*L + S);
    case FF_DIPOLE:
        return std::pow((L*L + mN*mN) / (L*L + S), 2.0);
    case FF_EXP_S:
        return std::exp((mN*mN - S) / (L*L));
    case FF_COV_DIPOLE:
        return std::pow((L*L - mMeson*mMeson) / (L*L - t), 2.0);
    case FF_DIPOLE_S_CH:
        return (std::pow(L,4) + std::pow(mBaryon,4)) /
               (std::pow(L,4) + sM*sM);
    case FF_EXP_T:
        return std::exp((t - mMeson*mMeson) / (L*L));
    case FF_PAULI_VILLAR: {
        double num = (t - mMeson*mMeson);
        double den = (t - L*L);
        return std::sqrt(std::max(0.0, 1.0 - (num*num)/(den*den)));
    }
    default:
        throw std::invalid_argument("Unknown FormFactorType");
    }
}

// ----------------------------------------------------------
//  fypiN  -- pi-N splitting function
//  Translated from FUNCTION fypiN(y,kT,L,typ,dis)
// ----------------------------------------------------------
inline double fypiN(double y, double kT, double L,
                    FormFactorType typ, DisChannel dis)
{
    if (y <= 0.0 || y >= 1.0) return 0.0;

    double g_piNN, gg;
    // dis=0: charge exchange  N -> pi^- + p
    // dis=1: neutral          N -> pi^0 + N
    if (dis == DIS_CHARGE_EXCHANGE) {
        g_piNN = std::sqrt(14.40 * 4.0 * PI);
        gg = 2.0 * g_piNN*g_piNN / (16.0 * PI*PI);
    } else {
        g_piNN = std::sqrt(14.40 * 4.0 * PI);
        gg =       g_piNN*g_piNN / (16.0 * PI*PI);
    }

    const double mP  = mN;   // intermediate baryon mass = nucleon
    const double kT2 = kT*kT;
    const double SpiN = (kT2 + mpi*mpi)/y + (kT2 + mP*mP)/(1.0-y);

    double t  = 0.0, sM = 0.0;
    if (typ == FF_COV_DIPOLE || typ == FF_EXP_T || typ == FF_PAULI_VILLAR)
        t = (-kT2 - mN*mN*y*y) / (1.0-y);
    if (typ == FF_DIPOLE_S_CH)
        sM = (kT2 + (1.0+y)*mpi*mpi)/y + (kT2 + y*mP*mP)/(1.0-y) + mN*mN;

    double FF = formFactor(SpiN, t, mpi, mP, L, typ, sM);

    double ss = (kT2 + std::pow(mP - (1.0-y)*mN, 2.0)) / (1.0-y)
              / std::pow((1.0-y)*(SpiN - mN*mN), 2.0)
              * FF*FF;

    return gg * (1.0-y) / y * ss;
}

// ----------------------------------------------------------
//  f_rhoN  -- rho-N splitting function
//  Translated from FUNCTION f_rhoN(y,kT,L,typ,dis)
// ----------------------------------------------------------
inline double f_rhoN(double y, double kT, double L,
                     FormFactorType typ, DisChannel dis)
{
    if (y <= 0.0 || y >= 1.0) return 0.0;

    double g_RoNN, f_RoNN, gg, fff, fg;
    // dis=0: N -> rho^- + proton  (factor 2 in g)
    // dis=1: N -> rho^0 + N       (factor 1 in g)
    if (dis == DIS_CHARGE_EXCHANGE) {
        g_RoNN = std::sqrt(2.0 * 0.55 * 4.0 * PI);
    } else {
        g_RoNN = std::sqrt(      0.55 * 4.0 * PI);
    }
    f_RoNN = 6.1 * g_RoNN;
    gg  = g_RoNN*g_RoNN / (16.0 * PI*PI);
    fff = f_RoNN*f_RoNN / (16.0 * PI*PI);
    fg  = f_RoNN*g_RoNN / (16.0 * PI*PI);

    const double mP   = mN;
    const double kT2  = kT*kT;
    const double SRoN = (kT2 + mrho*mrho)/y + (kT2 + mP*mP)/(1.0-y);

    double t = 0.0, sM = 0.0;
    if (typ == FF_COV_DIPOLE)
        t = (-kT2 - mN*mN*y*y) / (1.0-y);
    if (typ == FF_DIPOLE_S_CH)
        sM = (kT2 + (1.0+y)*mrho*mrho)/y + (kT2 + y*mP*mP)/(1.0-y) + mN*mN;

    double FF = formFactor(SRoN, t, mrho, mP, L, typ, sM);

    // TOPT with P[alpha]-p[alpha] derivative coupling
    double P_k  = (mrho*mrho + y*y*mN*mN + kT2) / (2.0*y);
    double P_p  = (mP*mP + (1.0-y)*(1.0-y)*mN*mN + kT2) / (2.0*(1.0-y));
    double pl_k = (mP*mP + kT2)*y / (2.0*(1.0-y))
                + (mrho*mrho + kT2)*(1.0-y) / (2.0*y)
                + kT2;

    double sv = -6.0*mN*mP + 4.0*P_k*pl_k/(mrho*mrho) + 2.0*P_p;

    double st = -(P_p*P_p)
              + P_p*(mP+mN)*(mP+mN)
              - mP*mN*(mP*mP + mN*mN + mP*mN)
              + 1.0/(2.0*mrho*mrho)
              * ( (P_p - mP*mN)*(P_k-pl_k)*(P_k-pl_k)
                - 2.0*(P_k-pl_k)*(mP*mP*P_k - mN*mN*pl_k)
                + 2.0*P_k*pl_k*(2.0*P_p - mP*mP - mN*mN) );

    double si = -4.0*(mP+mN)*(mP*mN - P_p)
              - 2.0*(mP*P_k*P_k - (mP+mN)*P_k*pl_k + mN*pl_k*pl_k)/(mrho*mrho);

    double ss = (gg*sv + fff*st/(mN*mN) + fg*si/mN)
              / std::pow(y*(SRoN - mN*mN), 2.0)
              * FF*FF;

    return y / (1.0-y) * ss;
}

// ----------------------------------------------------------
//  f_RhoDel  -- rho-Delta splitting function
//  Translated from FUNCTION f_RhoDel(y,kT,L,typ,dis)
// ----------------------------------------------------------
inline double f_RhoDel(double y, double kT, double L,
                       FormFactorType typ, DisChannel dis)
{
    if (y <= 0.0 || y >= 0.999) return 0.0;

    double g_NDS, gg;
    if (dis == DIS_CHARGE_EXCHANGE) {
        g_NDS = std::sqrt(20.448 * 4.0 * PI);
        gg = (2.0/3.0) * g_NDS*g_NDS / (16.0 * PI*PI * mrho*mrho);
    } else {
        g_NDS = 0.0;
        gg = g_NDS*g_NDS / (16.0 * PI*PI * mrho*mrho);
    }

    const double kT2  = kT*kT;
    const double SRoD = (kT2 + mrho*mrho)/y + (kT2 + mDelta*mDelta)/(1.0-y);

    double t = 0.0, sM = 0.0;
    if (typ == FF_COV_DIPOLE)
        t = (-kT2 - (1.0-y)*(mrho*mrho - y*mN*mN)) / y;
    if (typ == FF_DIPOLE_S_CH)
        sM = (kT2 + (2.0-y)*mrho*mrho)/(1.0-y)
           + (kT2 + (1.0-y)*mDelta*mDelta)/y
           + mN*mN;

    double FF = formFactor(SRoD, t, mrho, mDelta, L, typ, sM);

    double P_k  = (mrho*mrho + y*y*mN*mN + kT2) / (2.0*y);
    double P_p  = (mDelta*mDelta + (1.0-y)*(1.0-y)*mN*mN + kT2) / (2.0*(1.0-y));
    double pl_k = (mDelta*mDelta + kT2)*y / (2.0*(1.0-y))
                + (mrho*mrho + kT2)*(1.0-y) / (2.0*y)
                + kT2;

    double sr =
        -4.0*mN*mDelta/3.0*(2.0*mDelta*mDelta + mN*mDelta + 2.0*mN*mN)
        -4.0*mN*mDelta/(3.0*mrho*mrho)*std::pow(P_k - pl_k, 2.0)
        -4.0/(3.0*mrho*mrho)*(mDelta*mDelta*P_k*P_k + mN*mN*pl_k*pl_k)
        +4.0*P_p/3.0*(2.0*mDelta*mDelta + 4.0*mN*mDelta + mN*mN)
        +4.0*P_p/(3.0*mrho*mrho)*pl_k*pl_k*(1.0 - mN*mN/(mDelta*mDelta))
        -4.0*P_p*P_p*(1.0
                      - 2.0*P_k*pl_k/(3.0*mrho*mrho*mDelta*mDelta)
                      - P_p/(3.0*mDelta*mDelta));

    double ss = sr / std::pow((1.0-y)*(SRoD - mN*mN), 2.0) * FF*FF;

    return gg * (1.0-y) / y * ss;
}

// ----------------------------------------------------------
//  fypiD  -- pi-Delta splitting function
//  Translated from FUNCTION fypiD(y,kT,L,typ,dis)
// ----------------------------------------------------------
inline double fypiD(double y, double kT, double L,
                    FormFactorType typ, DisChannel dis)
{
    if (y <= 0.0 || y >= 0.999) return 0.0;

    double g_DLcN, gg;
    if (dis == DIS_CHARGE_EXCHANGE) {
        g_DLcN = std::sqrt(0.2237 * 4.0 * PI);
        gg = (2.0/3.0) * g_DLcN*g_DLcN / (16.0 * PI*PI);
    } else {
        g_DLcN = 0.0;
        gg = g_DLcN*g_DLcN / (16.0 * PI*PI);
    }

    const double kT2  = kT*kT;
    const double SpiD = (kT2 + mpi*mpi)/y + (kT2 + mDelta*mDelta)/(1.0-y);

    double t = 0.0, sM = 0.0;
    if (typ == FF_COV_DIPOLE)
        t = (-kT2 - mN*mN*y*y) / (1.0-y);
    if (typ == FF_DIPOLE_S_CH)
        sM = (kT2 + (1.0+y)*mpi*mpi)/y
           + (kT2 + y*mDelta*mDelta)/(1.0-y)
           + mN*mN;

    double FF = formFactor(SpiD, t, mpi, mDelta, L, typ, sM);

    double yR = 1.0 - y;
    double ss = (kT2 + std::pow(mDelta - yR*mN, 2.0))
              * std::pow(kT2 + std::pow(mDelta + yR*mN, 2.0), 2.0)
              / (6.0*mDelta*mDelta*std::pow(yR, 3.0))
              / std::pow((1.0-yR)*(SpiD - mN*mN), 2.0)
              * FF*FF;

    return gg / (mpi*mpi) * (1.0-yR) / yR * ss;
}

} // namespace PionCloud
