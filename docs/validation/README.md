# Validation against the published TDIS figures

Reproduces the pion-cloud figures of the JLab TDIS proposal and overlays PionCloud
on the original curves. The output is the figures in `docs/figures/` used in
Sec. "Worked examples" of the manual.

| Figure (PR12-15-006) | Mode  | Content |
|----------------------|-------|---------|
| Fig. 47 (left/right) | theta | F2(θp′) for p→π⁰p, ρ⁰p; Δx=[0,0.6], Δk=[0,500] MeV / Δx=[0.05,0.6], Δk=[60,250] MeV |
| Fig. 46              | x     | F2(x) for p→π⁰p, ρ⁰p; Δk=[0,500] MeV, all angles |
| Fig. 5 (= Hobbs Fig. 2R) | kbin | p→π⁰p in |k| bins 60–100, 100–200, 200–300, 300–400 MeV; total F2^(πp), F2p |
| Fig. 6               | kbin  | n→π⁻p, ρ⁻p, πΔ; Δk=[60,250] MeV, Δθ=[30°,180°]; total F2^(πN), F2n |
| Fig. 37 (= Hobbs Fig. 2L) | kbin | same channels; Δk=[250,400] MeV, Δθ=[30°,70°] |

## Sources

- A. Camsonne et al., JLab proposal PR12-15-006 (PAC43, 2015; approved as E12-15-006),
  https://www.jlab.org/exp_prog/proposals/15/PR12-15-006.pdf (pages 14, 15, 54, 66, 67).
- T. J. Hobbs, Few-Body Syst. 56, 363 (2015), arXiv:1411.2216, Fig. 2. Its curves are
  identical to Figs. 5 and 37 above.

## Files

- `ref_data/figN_<panel>_<colour>[_dash].dat`: columns `x  F2` (or `theta_deg  F2`). Each
  file is one published curve, extracted exactly from the vector graphics of the proposal
  PDF. The extraction converts the page with `pdftocairo -svg` and maps it with the axis
  frame and the major ticks, using `extract_curves.py`. No hand digitisation was involved.
  Colours follow the published figure.
- `extract_curves.py`: the extractor. Example:
  `pdftocairo -svg -f 67 -l 67 PR12-15-006.pdf p67.svg;
   extract_curves.py p67.svg fig47 200,330:left:0:15:0:-6:1:1 430,330:right:0:15:0:-7:1:1`
- `make_runs.sh [binary] [workdir]`: writes one input file per curve into
  `<workdir>/inputs/`, then runs PionCloud in `<workdir>/runs/`.
- `compare.py [workdir]`: makes the overlays (PyROOT) in `<workdir>/plots/` and the ratios
  in `<workdir>/summary.txt`.
- `run_compare.sh [binary] [workdir]`: does both.

```bash
source /path/to/root/bin/thisroot.sh       # PyROOT + numpy
docs/validation/run_compare.sh build/PionCloud /tmp/pc_validation
cp /tmp/pc_validation/plots/{fig47,fig46,fig6,fig37}.pdf docs/figures/
cp /tmp/pc_validation/plots/fig5_th30-70.pdf docs/figures/fig5.pdf
```

## Result (PionCloud v1.3, 2026-11-12)

- Tagged curves: median deviation +0.2% to +1.7%.
- Bulk (curve within ×100 of its maximum): at most 2–4%, and 7% for Fig. 6 just before
  the kinematic end point.
- Inclusive curves agree to 0.1% once the sum over charge states is taken:
  ×3 for neutral, ×3/2 for charge exchange.
- CTEQ6M F2p/F2n lies 3–7% above the published orange curves (open item).
- θe = 35° reproduces the figures. θe = 12° does not: it is 3–11% high at x ≳ 0.15.

## Automatic test (v1.4)

`ctest` runs `PionCloudValidation` (`tests/validation.cc`), which recomputes these curves and
compares them with `ref_data/`. Pass criteria:
- tagged curves: −1% < median < +3% and max |dev| < 8% in the bulk;
- total F₂^(πN) (`charge_states = total`): |median| < 0.5%.

With `kbin_method = exact` the kbin maxima halve: Fig. 6 6–7% → 3.5%, Fig. 37 ~2.5% → 1–1.6%.
