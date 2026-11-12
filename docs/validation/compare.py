#!/usr/bin/env python3
"""Overlay PionCloud output on the curves extracted from the published figures
(ref_data/, see README.md) and print PionCloud/reference ratios.

Usage: python3 compare.py [work dir]      (work dir as given to make_runs.sh)
Writes <work dir>/plots/<fig>.pdf/.png and <work dir>/summary.txt."""
import glob, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

WORK = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else 'validation_runs')
R, P = f'{WORK}/runs', f'{WORK}/plots'
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ref_data')
os.makedirs(P, exist_ok=True)
COL = {'black': ROOT.kBlack, 'red': ROOT.kRed+1, 'blue': ROOT.kBlue+1, 'green': ROOT.kGreen+2,
       'orange': ROOT.kOrange-3, 'violet': ROOT.kViolet+1}
KEEP = []   # keep ROOT objects alive

def latex(t):
    return (t.replace('$', '').replace('\\', '#').replace('Δ', '#Delta').replace('θ', '#theta')
             .replace('#to', '#rightarrow').replace('°', '#circ')
             .replace('^0', '^{0}').replace('^-', '^{-}').replace('F_2^n', 'F_{2}^{n}')
             .replace('F_2^p', 'F_{2}^{p}').replace('F_2^{', 'F_{2}^{'))

def ref(name):
    a = np.loadtxt(f'{D}/{name}.dat')
    return a[np.argsort(a[:, 0])]

def run(name, col=2, scale=1.0):
    f = glob.glob(f'{R}/{name}_xscan.dat') or glob.glob(f'{R}/{name}_thetascan.dat')
    a = np.loadtxt(f[0], comments='#')
    return a[:, 0], scale*a[:, col]

summary = []

def ratio_stats(label, xr, yr, xp, yp):
    """PionCloud/reference evaluated at the reference points (log-interpolated)."""
    m = (yp > 0)
    xp, yp = xp[m], yp[m]
    sel = (xr >= xp.min()) & (xr <= xp.max()) & (yr > 0)
    if sel.sum() < 2:
        return None
    yi = np.exp(np.interp(xr[sel], xp, np.log(yp)))
    r = yi / yr[sel]
    summary.append(f'  {label:42s} PC/ref: median {np.median(r):.3f}  range [{r.min():.3f}, {r.max():.3f}]  ({sel.sum()} pts)')
    return xr[sel], r

def graph(x, y, colour, dash, width, alpha=1.0):
    x = np.ascontiguousarray(x, dtype=float); y = np.ascontiguousarray(y, dtype=float)
    g = ROOT.TGraph(len(x), x, y)
    g.SetLineColorAlpha(COL[colour], alpha); g.SetMarkerColor(COL[colour])
    g.SetLineWidth(width); g.SetLineStyle(2 if dash else 1)
    KEEP.append(g)
    return g

def panel(c, i, n, curves, xlabel, title, ylim, xlim, legpos=(0.52, 0.52, 0.95, 0.88)):
    """curves: list of (label, refname, (runname, col), colour, dashed)"""
    x0, x1 = i/n, (i+1)/n
    pt = ROOT.TPad(f'pt{i}', '', x0, 0.30, x1, 1.0); pb = ROOT.TPad(f'pb{i}', '', x0, 0.0, x1, 0.30)
    for p in (pt, pb):
        p.SetLeftMargin(0.15); p.SetRightMargin(0.04); c.cd(); p.Draw(); KEEP.append(p)
    pt.SetBottomMargin(0.015); pt.SetTopMargin(0.08); pb.SetTopMargin(0.02); pb.SetBottomMargin(0.32)
    pt.cd(); pt.SetLogy(); pt.SetTicks(1, 1)
    fr = pt.DrawFrame(xlim[0], ylim[0], xlim[1], ylim[1]); fr.GetYaxis().SetTitle('F_{2}')
    fr.GetXaxis().SetLabelSize(0); fr.GetYaxis().SetTitleSize(0.055); fr.GetYaxis().SetLabelSize(0.045)
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.15, 0.94, latex(title)); KEEP.append(t)
    leg = ROOT.TLegend(*legpos); leg.SetBorderSize(0); leg.SetFillColorAlpha(ROOT.kWhite, 0.85); leg.SetTextSize(0.040)
    leg.SetHeader('thick: published / thin: PionCloud', ''); KEEP.append(leg)
    pb.cd(); pb.SetTicks(1, 1); pb.SetGridy()
    fb = pb.DrawFrame(xlim[0], 0.8, xlim[1], 1.2); KEEP.append(fb)
    fb.GetXaxis().SetTitle(latex(xlabel)); fb.GetYaxis().SetTitle('PionCloud / ref')
    for ax, s_ in ((fb.GetXaxis(), 0.12), (fb.GetYaxis(), 0.10)):
        ax.SetTitleSize(s_); ax.SetLabelSize(0.10)
    fb.GetYaxis().SetTitleOffset(0.6); fb.GetYaxis().SetNdivisions(504); fb.GetXaxis().SetTitleOffset(1.1)
    for lab, rn, spec, c_, dash in curves:
        pn = spec[0]
        if rn:
            a = ref(rn); pt.cd(); graph(a[:, 0], a[:, 1], c_, dash, 7, 0.35).Draw('L')
        if pn:
            xp, yp = run(*spec); m = yp > 0
            g = graph(xp[m], yp[m], c_, dash, 2); pt.cd(); g.Draw('L'); leg.AddEntry(g, latex(lab), 'l')
            if rn:
                rr = ratio_stats(f'{title}: {lab}', a[:, 0], a[:, 1], xp, yp)
                if rr is not None:
                    gr = graph(rr[0], rr[1], c_, dash, 2); gr.SetMarkerStyle(20); gr.SetMarkerSize(0.5)
                    pb.cd(); gr.Draw('LP')
    pt.cd(); leg.Draw()

def figure(name, panels, width=6.0):
    n = len(panels)
    c = ROOT.TCanvas(name, name, int(560*n), 560); KEEP.append(c)
    for i, pa in enumerate(panels):
        panel(c, i, n, *pa)
    for ext in ('pdf', 'png'):
        c.SaveAs(f'{P}/{name}.{ext}')
    summary.append('')

# ---------------- Fig. 47 ----------------
summary.append('Fig. 47 (PR12-15-006): theta_p\' dependence, p -> pi0 p / rho0 p')
figure('fig47', [
    ([('$\\pi^0 p$', 'fig47_left_red', ('fig47L_pion', 1), 'red', False),
      ('$\\rho^0 p$', 'fig47_left_blue_dash', ('fig47L_rho', 1), 'blue', True)],
     r"$\theta_{p'}$ [deg]", 'Fig. 47 left: Δx=[0,0.6], Δk=[0,500] MeV', (1e-6, 3e-3), (0, 90), (0.48, 0.04, 0.95, 0.30)),
    ([('$\\pi^0 p$', 'fig47_right_red', ('fig47R_pion', 1), 'red', False),
      ('$\\rho^0 p$', 'fig47_right_blue_dash', ('fig47R_rho', 1), 'blue', True)],
     r"$\theta_{p'}$ [deg]", 'Fig. 47 right: Δx=[0.05,0.6], Δk=[60,250] MeV', (1e-7, 2e-4), (0, 90), (0.48, 0.04, 0.95, 0.30)),
], width=5.5)

# ---------------- Fig. 46 ----------------
summary.append('Fig. 46 (PR12-15-006): x dependence, Δk=[0,500] MeV, all angles')
figure('fig46', [
    ([('$\\pi^0 p$', 'fig46_red', ('fig46_pion', 2), 'red', False),
      ('$\\rho^0 p$', 'fig46_blue_dash', ('fig46_rho', 2), 'blue', True)],
     '$x$', 'Fig. 46: Δk=[0,500] MeV, Δθ=[0,180]°', (1e-5, 1e-2), (0, 0.4), (0.15, 0.04, 0.55, 0.30)),
])

# ---------------- Fig. 5 / Hobbs Fig. 2 right ----------------
for th in ('30-180', '30-70'):
    summary.append(f'Fig. 5 (PR12-15-006) = Hobbs Fig. 2R: p -> pi0 p, |k| bins, angular window {th} deg')
    figure(f'fig5_th{th}', [
        ([('[60,100] MeV', 'fig5_main_black', (f'fig5_k060_th{th}', 2), 'black', False),
          ('[100,200] MeV', 'fig5_main_red_dash', (f'fig5_k100_th{th}', 2), 'red', True),
          ('[200,300] MeV', 'fig5_main_green_dash', (f'fig5_k200_th{th}', 2), 'green', True),
          ('[300,400] MeV', 'fig5_main_blue', (f'fig5_k300_th{th}', 2), 'blue', False),
          ('inclusive $F_2^{(\\pi p)}$ #times 3', 'fig5_main_violet_dash', ('fig5_inclusive', 2, 3.0), 'violet', True),
          ('$F_2^p$ (CTEQ6M)', 'fig5_main_orange', ('fig5_inclusive', 4), 'orange', False)],
         '$x$', f'Fig. 5: $p\\to\\pi^0 p$, Δθ=[{th.replace("-", ",")}]°', (1e-6, 1), (0.04, 0.3), (0.55, 0.40, 0.95, 0.78)),
    ])

# ---------------- Fig. 6 ----------------
summary.append('Fig. 6 (PR12-15-006): n -> pi- p etc., Δθ=[30,180] deg, Δk=[60,250] MeV assumed')
figure('fig6', [
    ([('$\\pi N$', 'fig6_main_black', ('fig6_pion', 2), 'black', False),
      ('$\\rho N$', 'fig6_main_red_dash', ('fig6_rho', 2), 'red', True),
      ('$\\pi\\Delta$', 'fig6_main_green_dash', ('fig6_pi_delta', 2), 'green', True),
      ('inclusive $F_2^{(\\pi N)}$ #times 3/2', 'fig6_main_violet_dash', ('fig6_inclusive', 2, 1.5), 'violet', True),
      ('$F_2^n$ (CTEQ6M)', 'fig6_main_orange', ('fig6_inclusive', 4), 'orange', False)],
     '$x$', 'Fig. 6: $n\\to\\pi^- p$, Δθ=[30,180]°', (5e-7, 1), (0.04, 0.2), (0.55, 0.40, 0.95, 0.78)),
])

# ---------------- Fig. 37 / Hobbs Fig. 2 left ----------------
summary.append('Fig. 37 (PR12-15-006) = Hobbs Fig. 2L: n -> pi- p etc., Δk=[250,400] MeV, Δθ=[30,70] deg')
figure('fig37', [
    ([('$\\pi N$', 'fig37_main_black', ('fig37_pion', 2), 'black', False),
      ('$\\rho N$', 'fig37_main_red_dash', ('fig37_rho', 2), 'red', True),
      ('$\\pi\\Delta$', 'fig37_main_green_dash', ('fig37_pi_delta', 2), 'green', True),
      ('inclusive $F_2^{(\\pi N)}$ #times 3/2', 'fig37_main_violet_dash', ('fig37_inclusive', 2, 1.5), 'violet', True),
      ('$F_2^n$ (CTEQ6M)', 'fig37_main_orange', ('fig37_inclusive', 4), 'orange', False)],
     '$x$', 'Fig. 37: $n\\to\\pi^- p$, Δk=[250,400] MeV, Δθ=[30,70]°', (5e-7, 1), (0.04, 0.3), (0.55, 0.40, 0.95, 0.78)),
])

open(f'{WORK}/summary.txt', 'w').write('\n'.join(summary) + '\n')
print('\n'.join(summary))
