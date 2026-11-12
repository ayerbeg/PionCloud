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

### Default model: `exp_s` with Λ = 1.56 GeV (was Pauli–Villars with Λ = 1.33 GeV)

**What.** `RunParams` defaults are now `typ = FF_EXP_S`, `L = 1.56` (`Ld = 1.39`
unchanged), in `PhysicsParams.hh`.

**Why.**
- The v1.3 pairing, inherited from `3Var_x.f`, matches no published fit. In McKenney et al.,
  PRD 93 (2016) 054011, Table I, **1.33 GeV is the best fit of the s-dependent exponential**,
  while the Pauli–Villars fit is **0.27 GeV**. The comment next to `L = 1.33` in `3Var_x.f`
  itself says "s-dependent exp. – PRD93,054011".
- The default `Ld = 1.39` was already the HSS `exp_s` value, so the defaults mixed two models.
- The C interface default (`pioncloud_init_default`) already used `exp_s` / 1.56. The C++ and
  C/Python defaults gave different physics.
- HSS `exp_s` (Holtmann–Szczurek–Speth, Λ_πN = Λ_ρN = 1.56 ± 0.07 GeV,
  Λ_πΔ = Λ_ρΔ = 1.39 ± 0.07 GeV, as in PR12-15-006) is the model that reproduces all the
  published TDIS figures (manual, Sec. 7).
- The PV implementation itself is correct: it is the effective form factor, Eq. (9) of
  McKenney et al.

**Numbers.**
- Change **only for runs that did not set `form_factor` and `lambda`**. On the Fig. 37
  configuration the old default gave F₂^(πN) **1.63–1.89×** the HSS model. Results produced
  with defaults should be checked.
- Runs with explicit settings: bit-identical (23/23).

**Also.** The constructor prints a warning when the cut-off used by the selected channel is far
outside the published range of its regulator. Ranges: PV 0.15–0.6, `exp_s` 0.9–2.0,
`exp_t` 0.5–1.5 GeV. The old default triggers it.

### Odd numbers of integration steps are now integrated correctly

**What.** `simpsonWeight(i, n)` in `Kinematics.hh` works for any n. For even n it is the
standard composite Simpson rule, unchanged. For odd n it uses Simpson on the first n−3
intervals plus Simpson's 3/8 rule on the last three (also exact for cubics). n = 1 is the
trapezoid rule.

**Why.** Simpson's 1-4-2-…-4-1 pattern is only valid for even n. v1.3 silently accepted odd
`nk`, `nx_int`, `ny`, `nkT`, … and returned a slightly wrong integral (θ scan with `nk = 201`:
−0.4%).

**Numbers.**
- Even n (all defaults and all shipped inputs): bit-identical.
- Odd n: now agrees with even n (`nk = 201` vs `200`: 0.000%).
- Tests check ∫1 and ∫x³ for n = 1, 3, 5, 7, 101.

### Δ channels with `dissociation = neutral` are rejected instead of returning 0

**What.** `InputReader::validate` and the `Calculator` constructor reject `pi_delta` or
`rho_delta` together with `dissociation = neutral`.

**Why.** In the Fortran (and in `fypiD`, `f_RhoDel`) the coupling for that combination is
zero. The run "succeeded" and silently produced F₂ = 0. For the Δ channels, `charge_exchange`
means the Fortran's sum over charge states (isospin factor 2/3, e.g. π⁰Δ⁰ + π⁻Δ⁺ for the
neutron, PR12-15-006 Fig. 6). This is documented at `DisChannel` in `PhysicsParams.hh`.

**Numbers.** No valid configuration changes.

### ρΔ channel selectable (`contribution = rho_delta`, `--flag 3`)

**What.** New `FLAG_RHO_DELTA = 3`, dispatched to the existing `f_RhoDel()`, which uses
`lambda_delta`. Added to input, CLI, output names (`rhoDelta`) and the C API.

**Why.** The Fortran ρΔ splitting function was ported in v1.3 but could not be selected.

**Caveat.** Not validated against a published tagged curve; none exists in the TDIS proposal.
On the Fig. 37 configuration it is 4–13% of πΔ, rising with x.

### Strange quarks in F₂^π (optional)

`pion_strange = true` in the input file (or `--strange`) adds 2/9·x s̄ to F₂^π. GRV92 has an
SU(3)-symmetric sea (s̄ = S); GRS99 has its own s̄.

**Why optional.** The Fortran omits the term (F₂^π = 5/9(xV + 2xS) is u, d only), and the
published figures are reproduced without it. Default `false`, so the default is
bit-identical.

### Command-line options that were documented but did not exist

`--grv`, `--nucleon`, `--precision`, `--inclusive`, `--verbose`, `--xmin-th`, `--xmax-th`,
`--kmin-th` and `--kmax-th` were listed by `--help` in v1.3, but every one of them exited with
"Unknown option". `TODO.md` listed them as done because their fields existed in `main.cc`
and were never parsed. They are implemented now, together with `--quiet` and `--strange`.

### Test suite: v1.3's ctest "Passed" did not mean the tests ran

**What was wrong.**
- In builds with CTEQ6, `PionCloudTests` ran in `build/tests/`, where the CTEQ6 tables are not
  present.
- At the first `Calculator` the Fortran `SetCtq6` could not open `cteq6m.tbl` and executed
  `STOP`, which ends the process with **exit status 0**.
- ctest therefore reported "Passed" after about 80 of the tests. The Calculator, kbin-scan,
  Fortran-regression, output and CTEQ6 tests never ran.

**Fix.**
1. `tests/CMakeLists.txt` runs the suite in the build directory, where the tables are, and
   accepts it only if the summary line `Results: N passed, 0 failed` appears. Any early exit
   or crash is now a failure.
2. The `Calculator` constructor checks that `cteq6m.tbl` is readable **before** calling the
   Fortran, and throws `std::runtime_error` with instructions instead of letting `STOP`
   terminate the host program. That matters for generators and Python embedding the library.
   `PionCloud` (main) reports it as `Fatal: …`.

**Result.** All **197** tests now run and pass, against about 80 actually executed before.
A missing table now gives exit status 134 (test binary) or 1 (`PionCloud`), and the ctest
regex no longer matches.

### Comment fix

`CTEQ6::xfx()` returns f(x, Q), not x·f. `Ctq6Pdf` "returns the parton distribution". The
F₂ formulas were already correct (they multiply by x); only the comment was wrong.

### Pluggable models (`include/Models.hh`), for theorists' new models

**What.**
- New interfaces `SplittingModel` (f(y, k_T)) and `PionSFModel` (F₂^M(x, Q²)).
- Built-ins `BuiltinSplitting` (the Fortran ports, selected by `contribution`) and
  `BuiltinPionSF` (GRV92/GRS99).
- A `ModelRegistry` that selects models by name. Input keys: `splitting_model`,
  `pion_sf_model`.
- A second `Calculator` constructor that injects user models.

**Why.**
- The physics factorises as (splitting function) × (meson structure function) over a
  kinematic region.
- In v1.3 the first two were hard-wired with `switch` statements in four places. A new
  model therefore meant editing the integration loops.
- Now a theorist adds a model as one class, with no change to PionCloud's sources, and all
  scans, outputs and the C interface use it.
- Examples in the header of `Models.hh`; tests in `test_Models`. The tests check that an
  injected 2× model gives exactly 2× the result, that selection by name works, and that
  unknown names are rejected.

**Numbers.** Bit-identical (23/23). No measurable speed change.

### Code documentation

**What.** Physics block comments added:
- `SplittingFunctions.hh`: definition and normalisation of f_MB (⟨n⟩ = ∫dy ∫dk_T² f),
  variables, form-factor conventions (`exp_s` vs the HSS 2Λ² convention), couplings, isospin
  meaning of `dis` per channel, and "how to add a channel".
- `Kinematics.hh`: one angle, three names (θ_h = θ_p′ = φ_k); the y = z̃ formula with its
  source (Hobbs 2015, Eq. 5); the Q² formula.
- `Calculator.hh`: "how to modify"; which `RunParams` fields each scan reads.
- `GRV.hh`: valence/sea normalisation per set.
- New `CONTRIBUTING.md`: workflow, numerical-change policy, tests, adding physics.

**Verification of the documented normalisation** (independent publication). Summed over
charge states, PionCloud gives:
- `exp_s`, Λ = 1.33 GeV: ⟨n⟩_πN = **0.224**;
- Pauli–Villars, Λ = 0.27 GeV: ⟨n⟩_πN = **0.310**.

McKenney et al., PRD 93 054011, Table I, give 0.23 and 0.31.

**Open finding (not changed).** ⟨n⟩_πΔ = 0.027 at `exp_s` 1.33 GeV, against 0.06 in the same
table: a factor ≈ 2.2. PionCloud's πΔ reproduces Hobbs' Fortran and PR12-15-006 Fig. 6 to 1%,
so the difference lies between the Fortran πΔ and the 2016 paper. Candidates:
- the isospin factor 2/3;
- the coupling normalisation 0.2237;
- off-shell / δ-function terms counted in the paper.

Needs a term-by-term comparison with Eq. (5) of the paper (TODO).

**Numbers.** Comments only: bit-identical.

### kbin scan: exact integration region (`kbin_method = exact`)

**What.**
- `kbinF2Exact()` (`Calculator.cc`) and `include/Quadrature.hh`.
- With u = k_T², A = 2m_N(1−y) and c = y(2−y)m_N², the f2_pi_sub kinematics give
  |k|²(u) = u + (u+c)²/A² and cos φ_k = (u+c)/(A|k|).
- Each cut becomes a quadratic inequality in u:
  - the |k| window is an interval;
  - cos φ_k ≤ c_max lies between two roots;
  - cos φ_k ≥ c_min lies outside two roots.
- So the allowed region at fixed y is at most two u-intervals, and no y above
  y_max(k_max) contributes.
- Both are integrated with Gauss–Legendre: y in `gl_ypanels` panels; u in pieces split at fixed
  breakpoints, `gl_nodes` nodes per piece.
- The old algorithm is kept unchanged as `kbin_method = grid`, the default. It now lives in
  `kbinF2Grid()`, shared by `runKBinScan()` and `eval()`, which previously had separate copies
  of the loops.

**Why.**
- The f2_pi_sub method puts k_T on a uniform grid over [0, 10] GeV and rejects points outside
  the windows with `continue`. The integrand is then a step function, so Simpson converges only
  slowly.
- For |k| ≤ 0.4 GeV more than 96% of the k_T points are rejected.
- The error is dominated by the y grid near the kinematic end point.

**Evidence** (`analyses/code_review/kbin_convergence`). Maximum deviation from the exact method
at 32×64, over the x points within ×1000 of the maximum:

| configuration | grid 100×5 000 (Fortran setting) | grid 400×200 000 (≈15 s) | exact 16×32 (0.3 s) |
|---|---|---|---|
| Fig. 5, \|k\| = 60–100 MeV | up to 100% (end point) | 7.8% | 0.0006% |
| Fig. 6 πN | 13.5% | 0.78% | 0.002% |
| Fig. 37 πN | 5.0% | 0.25% | 0.002% |

- The exact method agrees with itself at two resolutions to 0.002%. The grid method converges
  toward it.
- **Against the published curves**, the bulk maximum deviation halves with the exact method:
  Fig. 6 from 6–7% to 3.5%, Fig. 37 from about 2.5% to 1–1.6%. The largest deviations in the
  validation were therefore mostly grid artefacts.

**Tests.**
- Gauss–Legendre is exact up to degree 2n−1.
- The analytic region matches a brute-force cut test on 9,600 random points with 0 mismatches.
- Exact vs grid (400×50 000) agree to < 0.5%.
- `eval()` equals the scan.

**Numbers.** Default (`grid`): bit-identical (23/23). Recommendation: use `kbin_method = exact`
for new results. Keep `grid` to reproduce the Fortran.

### `q2_fixed` and `charge_states = total`

- **`q2_fixed = <GeV²>`.** A fixed Q² at every x, instead of Q²(x) at a fixed electron angle.
  Useful for comparisons with theory at a fixed scale. One function `Calculator::Q2at()` now
  serves all scans and `eval()`.
- **`charge_states = total`.** Returns the sum over charge states for πN and ρN: ×3 for neutral,
  ×3/2 for charge exchange (isospin 1 + 2). This is the "fully integrated" F₂^(πN) of the
  published figures, which v1.3 users had to scale by hand. The Δ channels are already summed
  in the Fortran, so their factor is 1.

**Tests.** Q² constant; total = 3× neutral = 3/2× charge exchange = neutral + charge exchange.
Defaults: bit-identical.

### Automatic validation against the published TDIS figures (ctest `PionCloudValidation`)

**What.** `tests/validation.cc` runs the 17 configurations of Figs. 5, 6, 37, 46 and 47
(PR12-15-006) and compares them with `docs/validation/ref_data`.

**Pass criteria.**
- Tagged curves: −1% < median < +3% and max |dev| < 8% in the bulk.
- Totals: |median| < 0.5%.

Measured: medians +0.3% to +1.7%, max ≤ 6.9%, totals −0.1%.

**Why.** It protects the agreement with the published physics against future changes; it runs
in 2 s. Shown to work: with Λ = 1.48 GeV (HSS −1σ), 15 of 17 checks fail. The two that pass
are πΔ, which does not use `lambda`.

### `kbin_method = exact` is now the DEFAULT; the original algorithm is selected with a flag

**Decision** (C. Ayerbe, 2026-11-12). Make the new integration the default, mark it clearly as
a new feature, and keep the original Fortran algorithm available behind an explicit switch.

**How to get the ORIGINAL Fortran (`f2_pi_sub.f90`) algorithm, i.e. v1.3 results:**

| interface | switch |
|-----------|--------|
| input file | `kbin_method = grid` |
| command line | `--kbin-grid` (or `--kbin-method grid`) |
| C / Python (ctypes) | `pioncloud_set_kbin_method(0)` (new function) |
| C++ | `params.kbin_method = PionCloud::KBIN_GRID;` |

**Where it is visible.**
- Every kbin run prints the method in the run summary ("kbin method: exact (NEW in v1.4, …)")
  and at the start of the scan.
- The method is also written into `<base>_params.txt` and the header of `<base>_xscan.dat`.

**Effect on results** (the 17 kbin validation runs, exact vs grid at the Fortran settings
ny = 100, nkT = 5000):
- median change per curve ≤ 0.3%;
- bulk (curve within ×100 of its maximum) typically 2–6%, near the end of the bulk;
- up to 35% for the narrow |k| = 60–100 MeV bin, where the grid result is not converged (see the
  convergence table above: up to 100% off at the Fortran settings).
- x- and θ-scans are not affected.
- With `kbin_method = grid`, all 23 validation runs are bit-identical to v1.3.

**Also affected.**
- `Calculator::eval()` and the C interface now use exact by default.
- `precision` presets and `ny`/`nkT` act only on the grid method. The exact method uses
  `gl_nodes` (16) and `gl_ypanels` (32).

**Tests.**
- The Fortran regression tests use `kbin_method = grid` explicitly, since their references are
  Fortran output. The same references are also checked with exact; both pass within 5%.
- `PionCloudValidation` checks every kbin figure with both methods: 28 checks. Exact: medians
  +0.3–1.7%, maxima ≤ 3.6%. Grid: maxima ≤ 6.9%.
- `docs/validation/make_runs.sh` sets `kbin_method = grid` explicitly, so the overlays in
  `docs/figures` stay those of the original algorithm.

### πΔ normalisation vs McKenney et al. (documented, not changed)

**Decision.** Keep the Fortran πΔ model and document the difference:
`docs/notes/piDelta_normalisation.pdf`, reproducible with `docs/notes/piDelta_check.cc`.

**Summary.** At `exp_s` with Λ = 1.33 GeV:
- PionCloud gives ⟨n⟩_πΔ = 0.027; McKenney et al. (PRD 93 054011, Table I) give 0.06.
- The paper's Eq. (5), implemented here with the same regulator, gives 0.058. The paper's
  number is therefore reproduced, and the difference is a difference of **model**, not a bug.

The factor 2.2 decomposes into:
- ×1.5 from isospin bookkeeping. This is the interpretation that the Fortran factor 2/3 is a
  partial charge-state sum; its exact meaning is undocumented.
- ×1.43 from the functional form: Fortran/HSS "Sullivan" Δ-pole form (∝ k_T⁶) vs the paper's
  chiral-EFT on-shell form (∝ k_⊥²).
- The couplings agree within 4%.
