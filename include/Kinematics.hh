#pragma once
// ============================================================
//  Kinematics.hh
//  Pion Cloud Model -- kinematic variable calculations
//
//  Provides the coordinate transforms between (k, theta_h)
//  and (y, kT) used in both 3Var_x and 3Var_theta modes,
//  as well as Q2 from beam kinematics.
// ============================================================

#include "PhysicsParams.hh"
#include <cmath>

namespace PionCloud {

// ----------------------------------------------------------
//  Q2 from inclusive electron kinematics
//  E     -- beam energy [GeV]
//  x     -- Bjorken-x
//  theta -- electron scattering angle [deg]
// ----------------------------------------------------------
inline double calcQ2(double E, double x, double theta_deg)
{
    // Guard: x=0 gives Q2=0 by definition; avoid division by zero.
    if (x <= 0.0) return 0.0;
    double th_rad = (theta_deg / 180.0) * PI;
    double sin2   = std::pow(std::sin(th_rad), 2.0);
    double denom  = (2.0*E / (x*mN)) * sin2 + 1.0;
    return 2.0 * mN * x * E * (1.0 - 1.0 / denom);
}

// ----------------------------------------------------------
//  (k, theta_h) --> (y, kT) coordinate transform
//  k      -- hadron momentum magnitude [GeV]
//  theta  -- hadron production angle   [rad]
//  y      -- IMF momentum fraction     [out]
//  kT     -- transverse momentum       [out] [GeV]
// ----------------------------------------------------------
inline void kThetaToYkT(double k, double theta,
                        double &y, double &kT)
{
    y  = (k/mN) * std::cos(theta)
       + (1.0/mN) * (mN - std::sqrt(mN*mN + k*k));
    kT = k * std::sin(theta);
}

// ----------------------------------------------------------
//  Phase-space factor (Jacobian) for the (k, theta_h) -> (y,kT)
//  transformation, as written in the Fortran programs.
//  rho_PS = (2/mN) * k^2 * sin(theta) * (1 - sin(phi_k)*cos(theta))
//  where phi_k = atan(k/mN)
// ----------------------------------------------------------
inline double phaseSpace(double k, double theta)
{
    double phi_k = std::atan(k / mN);
    return (2.0/mN) * k*k * std::sin(theta)
         * (1.0 - std::sin(phi_k)*std::cos(theta));
}

// ----------------------------------------------------------
//  Simpson's rule weight for step i (0-based)
//  Returns 1, 4, 2, 4, 2, ... (standard composite rule)
// ----------------------------------------------------------
inline double simpsonWeight(int i, int n)
{
    if (i == 0 || i == n) return 1.0;
    return (i % 2 == 1) ? 4.0 : 2.0;
}

} // namespace PionCloud
