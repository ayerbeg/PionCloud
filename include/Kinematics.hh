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
//  Composite quadrature weight for node i of n intervals (n+1 nodes,
//  i = 0..n), in units of h/3:   integral ~= (h/3) * sum_i w_i f_i
//
//  n even : standard composite Simpson     1, 4, 2, 4, ..., 2, 4, 1
//  n odd  : Simpson on [0, n-3] + Simpson's 3/8 rule on the last three
//           intervals, whose weights (3h/8)(1,3,3,1) are (h/3)(9/8,27/8,27/8,9/8).
//           Both pieces are exact for cubics, so the accuracy is the same
//           as for even n.
//  n == 1 : trapezoid, (h/2)(1,1) = (h/3)(3/2,3/2).
//
//  WHY (v1.4): the plain Simpson pattern is only valid for even n.  Up to
//  v1.3 an odd number of steps was accepted silently and gave a slightly
//  wrong integral (e.g. nk = 201 changed a theta scan by -0.4%).
//  For even n the weights are unchanged, so all previous results with
//  even step counts are reproduced bit for bit.
// ----------------------------------------------------------
inline double simpsonWeight(int i, int n)
{
    if (n % 2 == 0) {                       // ---- standard Simpson ----
        if (i == 0 || i == n) return 1.0;
        return (i % 2 == 1) ? 4.0 : 2.0;
    }
    if (n == 1) return 1.5;                 // ---- trapezoid ----
    const int m = n - 3;                    // Simpson part covers [0, m], m even
    double w = 0.0;
    if (i <= m) {                           // Simpson part
        if (i == 0 || i == m) w = 1.0;
        else                  w = (i % 2 == 1) ? 4.0 : 2.0;
        if (m == 0 && i == 0) w = 0.0;      // n == 3: no Simpson part
    }
    switch (i - m) {                        // 3/8 tail on [m, n]
    case 0: w += 9.0/8.0;  break;
    case 1: w += 27.0/8.0; break;
    case 2: w += 27.0/8.0; break;
    case 3: w += 9.0/8.0;  break;
    default: break;
    }
    return w;
}

} // namespace PionCloud
