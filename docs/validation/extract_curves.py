#!/usr/bin/env python3
"""Extract vector curves from the xmgrace figures of the TDIS proposal PR12-15-006
(page converted with `pdftocairo -svg`).

All paths are mapped to absolute page coordinates.  A panel is the largest black
rectangle (axis frame) containing a given page point.  The longest tick marks on the
frame are the major ticks; their values are supplied on the command line, which
calibrates the page -> data mapping (linear, or log10 for log axes).

Usage: extract_curves.py page.svg out_prefix spec [spec ...]
  spec = "px,py:name:x0:dx:xlog:y0:dy:ylog"
    px,py : any page point inside the panel frame (pt, SVG y downwards)
    x0,dx : value of the first (lowest) major x tick and the step between major ticks
            (for a log axis give log10 of the value and decades per major tick)
Writes <out_prefix>_<name>_<colour>[_dash][_n].dat, columns: x  y
"""
import re, sys, collections

COLS = [((1, 0, 0), 'red'), ((0, 0, 1), 'blue'), ((0, 0, 0), 'black'), ((0, .545, 0), 'green'),
        ((0, .392, 0), 'green'), ((1, .647, 0), 'orange'), ((.58, 0, .827), 'violet'),
        ((.627, .125, .941), 'violet'), ((1, 0, 1), 'magenta')]

def colname(st):
    v = [float(t)/100 for t in re.findall(r'([\d.]+)%', st)]
    return min(COLS, key=lambda c: sum((a-b)**2 for a, b in zip(c[0], v)))[1]

def load(svg):
    s = open(svg).read()
    out = []
    for p in re.findall(r'<path ([^>]*)/>', s):
        st = re.search(r'stroke="([^"]*)"', p)
        if not st:
            continue
        tr = re.search(r'transform="matrix\(([^)]*)\)"', p)
        a, b, c, d, e, f = map(float, tr.group(1).split(',')) if tr else (1, 0, 0, 1, 0, 0)
        P = [(float(x), float(y)) for x, y in re.findall(r'([-\d.e]+) ([-\d.e]+)', re.search(r' d="([^"]*)"', p).group(1))]
        P = [(a*x + c*y + e, b*x + d*y + f) for x, y in P]
        out.append((colname(st.group(1)), 'stroke-dasharray' in p, P))
    return out

def bb(P):
    xs = [x for x, _ in P]; ys = [y for _, y in P]
    return min(xs), max(xs), min(ys), max(ys)

def main():
    svg, out = sys.argv[1], sys.argv[2]
    allp = load(svg)
    for spec in sys.argv[3:]:
        pt, name, x0, dx, xlog, y0, dy, ylog = spec.split(':')
        px, py = map(float, pt.split(','))
        x0, dx, y0, dy = map(float, (x0, dx, y0, dy))
        frames = [bb(P) for c, _, P in allp if c == 'black' and len(P) in (4, 5)
                  and bb(P)[0] < px < bb(P)[1] and bb(P)[2] < py < bb(P)[3]]
        X0, X1, Y0, Y1 = min(frames, key=lambda t: (t[1]-t[0])*(t[3]-t[2]))
        tol = 0.02*max(X1-X0, Y1-Y0)
        segs = [bb(P) for c, _, P in allp if c == 'black' and len(P) == 2]
        # x ticks: vertical segments touching bottom frame edge; y ticks: horizontal touching left edge
        vt = [t for t in segs if t[1]-t[0] < 1e-3 and X0-1e-3 <= t[0] <= X1+1e-3 and abs(t[3]-Y1) < 1e-2 and t[3]-t[2] < 0.3*(Y1-Y0)]
        ht = [t for t in segs if t[3]-t[2] < 1e-3 and Y0-1e-3 <= t[2] <= Y1+1e-3 and abs(t[0]-X0) < 1e-2 and t[1]-t[0] < 0.3*(X1-X0)]
        Lv = max(round(t[3]-t[2], 2) for t in vt); Lh = max(round(t[1]-t[0], 2) for t in ht)
        xm = sorted(set(round(t[0], 2) for t in vt if abs(t[3]-t[2]-Lv) < 0.05))
        ym = sorted(set(round(t[2], 2) for t in ht if abs(t[1]-t[0]-Lh) < 0.05), reverse=True)  # page y down
        sx = (xm[-1]-xm[0])/(len(xm)-1); sy = (ym[0]-ym[-1])/(len(ym)-1)
        fx = lambda u: x0 + dx*(u-xm[0])/sx
        fy = lambda v: y0 + dy*(ym[0]-v)/sy
        print(f'{name}: frame x[{fx(X0):.4g},{fx(X1):.4g}] y[{fy(Y1):.4g},{fy(Y0):.4g}] '
              f'major ticks x={len(xm)} y={len(ym)}')
        n = collections.Counter()
        for c, dash, P in allp:
            if len(P) < 3:
                continue
            b = bb(P)
            if not (X0-tol <= b[0] and b[1] <= X1+tol and Y0-tol <= b[2] and b[3] <= Y1+tol):
                continue
            if c == 'black' and (abs(b[1]-b[0]) < 1e-2 or abs(b[3]-b[2]) < 1e-2 or len(P) <= 5 and b[1]-b[0] > 0.9*(X1-X0)):
                continue  # frame / ticks
            key = f'{c}{"_dash" if dash else ""}'
            n[key] += 1
            fn = f'{out}_{name}_{key}{"_"+str(n[key]) if n[key] > 1 else ""}.dat'
            with open(fn, 'w') as f:
                f.write(f'# digitised (vector) from {svg}, panel {name}, colour {c}, dashed={dash}\n')
                for u, v in P:
                    xv, yv = fx(u), fy(v)
                    f.write(f'{10**xv if xlog=="1" else xv:.6g} {10**yv if ylog=="1" else yv:.6g}\n')
            print(f'   {fn}: {len(P)} pts, x[{fx(b[0]):.3g},{fx(b[1]):.3g}]')

main()
