# PionCloud changelog

Each entry says **what** changed, **why**, **how it was verified**, and **whether numbers
change**. "Bit-identical" means that all 23 validation runs (`docs/validation/`, the inputs
for Figs. 5, 6, 37, 46 and 47 of the TDIS proposal) give exactly the same output as before
the change.

## v1.4 (in progress, branch `docs/manual-figures`)

### GRV99 fix: the pion PDF set selected by `grv_version = 99` was wrong

**What.** `include/GRV.hh`: `GRV99()` / `F2piGRV99()` now implement the LO parametrisation
of Glück, Reya, Schienbein, *Pionic parton distributions revisited*, Eur. Phys. J. C10 (1999)
313, arXiv:hep-ph/9903288, Appendix A ("GRS99"; the code keeps the historical name GRV99).
The new function `GRS99(x, Q2)` returns all four distributions in the paper's own
convention: total valence, light sea, strange sea and gluon.

**Why.** The v1.3 routine, a translation of `GRV99` in `f2_pi_sub.f90`, had three errors.
In the Fortran its call was commented out, so it was never used or validated there.

| # | v1.3 | GRS99 paper | effect |
|---|------|-------------|--------|
| 1 | valence `xV` used as **one** flavour | `v^π = u_v^π + d̄_v^π`, **both** valence quarks (text after Eq. 2; ∫v dx = 2) | valence part of F₂^π ×2; F₂^(πN) +107% on Fig. 37 |
| 2 | "sea" built from Eqs. (A.7)–(A.8) | those are the **strange** sea s^π; the light sea q̄^π is Eq. (A.4) with (A.6) | wrong sea shape/size |
| 3 | N = 1.212 + 0.498 s **−** 0.009 s² | N = 1.212 + 0.498 s **+** 0.009 s² (Eq. A.3) | ≲0.5% |

The per-flavour interface `GRV99(x, Q2, xV, xS)` now has the **same meaning as `GRV92`**:
xV is one valence flavour = v^π/2, and xS is one light-sea flavour = q̄^π.
F₂^π = 5/9 (xV + 2xS) is therefore correct for both sets.

**Evidence.**
- Sum rules computed from the code's own formulas (Q² = 1, 2, 4 GeV²):
  - GRV92: ∫V dx = 0.997 / 0.995 / 0.994.
  - v1.3 GRV99: ∫V dx = 1.996 / 1.989 / 1.983.
- After the fix, GRS99 satisfies the paper's sum rules: ∫v dx = 2.00, and the momentum sum
  ∫x(v + 4q̄ + 2s̄ + g) dx = 1.000 ± 0.005 for Q² = 1–100 GeV². If v were per flavour, the
  momentum sum would come out ≈ 1.4.
- The paper also gives ∫x v dx = 0.563 at μ², twice the per-flavour momentum.

**Tests.** `tests/test_all.cc::test_PionPDFSumRules` checks:
- the valence number for GRV92 and GRV99 (per flavour = 1);
- the GRS99 total valence (= 2) and momentum sum (= 1) at Q² = 1, 4 and 10 GeV²;
- the GRS99/GRV92 F₂^π ratio.

Run on the **v1.3** `GRV.hh`, the same assertions fail: ∫V = 1.986, F₂^π ratio = 1.84.

**Numbers.**
- Change **only** with `grv_version = 99`. On the Fig. 37 configuration (πN, Δk =
  250–400 MeV) F₂^(πN) goes from ≈ 2× GRV92 to 1.05–1.17× GRV92. That is a genuine
  PDF-set difference, largest at small x where the GRS99 sea differs from GRV92's zero
  input sea.
- Default runs (GRV92): bit-identical (23/23).
- Results produced with `grv_version = 99` before v1.4, including the 2020
  `TDISfunction/GRV91vsGRV99*.pdf` comparisons, should be redone.

**Also added.** An optional strange-quark term in F₂^π: `F2piGRV(x, Q2, true)` and
`F2piGRV99(x, Q2, true)` add 2/9·x s̄. It is off by default, which reproduces the Fortran.
It is exposed in the input file by the following entry.
