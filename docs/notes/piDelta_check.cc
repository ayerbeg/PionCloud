// Reproducer for docs/notes/piDelta_normalisation.tex (PionCloud v1.4).
// Build from the repository root:
//   g++ -std=c++17 -O2 -Iinclude docs/notes/piDelta_check.cc -o piDelta_check && ./piDelta_check
// (header-only: needs no library, no CTEQ6, no ROOT).  Runtime ~1 min.
// pi-Delta multiplicity: PionCloud/Fortran "Sullivan" form vs McKenney et al.
// PRD 93 054011 Eq. (5) (chiral-EFT on-shell form), same s-dependent exponential
// regulator F = exp[(M^2 - s)/Lambda^2], s = (kT^2+mpi^2)/y + (kT^2+MD^2)/(1-y).
#include "SplittingFunctions.hh"
#include "Kinematics.hh"
#include <cstdio>
using namespace PionCloud;
const double M = mN, MD = mDelta, mp = mpi;
// f_{pi+ Delta0}^{(on)}(y, kT) per dkT^2, Eq. (5); C = g^2/((4pi)^2 18 MD^2), g = 3 sqrt2/5 gA/fpi
double fMcK(double y, double kT, double L) {
    if (y <= 0 || y >= 1) return 0;
    const double gA = 1.267, fpi = 0.093, g = 3.0*std::sqrt(2.0)/5.0*gA/fpi;
    const double C = g*g / (16*PI*PI*18*MD*MD);
    const double k2 = kT*kT, Mb = MD + M, Dl = MD - M;
    const double D = -(k2 - y*(1-y)*M*M + y*MD*MD + (1-y)*mp*mp)/(1-y);   // Eq. (6)
    const double s = (k2 + mp*mp)/y + (k2 + MD*MD)/(1-y);
    const double F = std::exp((M*M - s)/(L*L));
    const double br = (Mb*Mb - mp*mp)*(Dl*Dl - mp*mp) - (3*(Dl*Dl - mp*mp) + 4*M*MD)*D;
    return C * y*(Mb*Mb - mp*mp)/((1-y)*D*D) * br * F*F;
}
template <class Fn> double I(Fn f) {
    const int ny = 2000, nu = 4000; const double umax = 25.0; double s = 0;
    for (int i = 0; i <= ny; ++i) { double y = 1e-6 + (1-2e-6)*i/ny, si = 0;
        for (int j = 0; j <= nu; ++j) { double u = umax*j/nu; si += simpsonWeight(j,nu)*f(y,std::sqrt(u)); }
        s += simpsonWeight(i,ny)*si*(umax/nu)/3; }
    return s*((1-2e-6)/ny)/3;
}
int main() {
    printf("g_piNDelta = %.2f GeV^-1 (paper: 11.8);  g^2/4pi = %.2f GeV^-2;  Fortran 0.2237/mpi^2 = %.2f GeV^-2\n",
           3*std::sqrt(2.)/5*1.267/0.093, std::pow(3*std::sqrt(2.)/5*1.267/0.093,2)/(4*PI), 0.2237/(mpi*mpi));
    for (double L : {1.33, 1.39}) {
        double nF = I([&](double y,double kT){ return fypiD(y,kT,L,FF_EXP_S,DIS_CHARGE_EXCHANGE); });
        double nM = 6*I([&](double y,double kT){ return fMcK(y,kT,L); });
        printf("exp_s Lambda = %.2f GeV:  PionCloud/Fortran (factor 2/3) <n> = %.4f ; same form x3/2 (all charges, total weight 1) = %.4f ; McKenney Eq.(5) x6 = %.4f\n",
               L, nF, 1.5*nF, nM);
    }
    // shape comparison at Lambda = 1.33: <y> and <kT^2>
    double L = 1.33;
    auto mom = [&](auto f){ double n=I(f), my=I([&](double y,double k){return y*f(y,k);}), mk=I([&](double y,double k){return k*k*f(y,k);}); printf("   <y> = %.3f  <kT^2> = %.3f GeV^2\n", my/n, mk/n); };
    printf("Fortran/Sullivan form:"); mom([&](double y,double k){return fypiD(y,k,L,FF_EXP_S,DIS_CHARGE_EXCHANGE);});
    printf("McKenney Eq.(5):      "); mom([&](double y,double k){return fMcK(y,k,L);});
}
