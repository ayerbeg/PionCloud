#pragma once
// ============================================================
//  SplittingFunctions.hh
//  Pion Cloud Model -- meson-baryon splitting functions
//
//  Direct C++ translations of the Fortran functions:
//    fypiN    -- pi-N splitting function
//    f_rhoN   -- rho-N splitting function
//    f_RhoDel -- rho-Delta splitting function
//    fypiD    -- pi-Delta splitting function
//
//  All functions are stateless and thread-safe.
//  Arguments:
//    y   -- IMF momentum fraction carried by the meson
//    kT  -- transverse momentum [GeV]
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
