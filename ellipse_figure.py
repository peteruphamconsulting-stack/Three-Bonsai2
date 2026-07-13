#!/usr/bin/env python3
"""
ellipse_figure.py -- render the certified least-prime-factor ELLIPSE packing (the prime cells
[p], areas d_p) as a flat figure, coloured by scale, with the largest cells labelled by prime.
Companion to ellipse_certify.py; reads the same head JSON and reconstructs each cell from its
stored centre, orientation, and area share.

Usage:
    python3 ellipse_figure.py --in ellipse_head_n1650.json --out ellipse_packing.png --dpi 300
    python3 ellipse_figure.py --in ellipse_head_n1650.json --label-n 40    # label the first 40 primes
    python3 ellipse_figure.py --in ellipse_head_n1650.json --label-n 0     # no labels
Requires: numpy, sympy, matplotlib.
"""
import argparse, json
from math import sqrt, log
import numpy as np, sympy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

RAMP = ["#0d2a33", "#164a52", "#1f6e6a", "#2f8f7d", "#57ab8e", "#8fc9a8"]   # deep -> light
EDGE = "#20302e"

def main():
    ap=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", required=True)
    ap.add_argument("--out", default="ellipse_packing.png"); ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--label-n", type=int, default=40,
                    help="label the cells of the first N primes (0 = no labels)")
    ap.add_argument("--font-base", type=float, default=6.0, help="min label font size (pt)")
    ap.add_argument("--font-scale", type=float, default=60.0, help="label font growth with cell size")
    ap.add_argument("--font-max", type=float, default=40.0, help="max label font size (pt)")
    a=ap.parse_args()
    d=json.load(open(a.infile)); V=d.get('verification', d)
    rho=float(V['aspect_ratio_b_over_a']); C=V['centers']
    ori=[1 if o==90 else 0 for o in V['orientation_pattern_degrees']]
    shares=V['area_shares_float']; n=len(shares)
    primes=list(sympy.primerange(2, sympy.prime(n)+1))[:n]         # cell i is the prime primes[i]
    sc=np.sqrt(np.array(shares))
    # colour on a log-log scale of the prime (like the rectangle figure): spreads the palette
    # across the whole packing instead of piling small cells into the light end.
    ll = np.log(np.log(np.array(primes, float)+1.0))
    lo, hi = ll.min(), ll.max()
    def col_i(i):
        t = (ll[i]-lo)/(hi-lo+1e-12)                    # 0 (p=2, deep) -> 1 (largest prime, light)
        return RAMP[min(len(RAMP)-1, int(t*len(RAMP)))]
    fig,ax=plt.subplots(figsize=(11, 11*rho*1.06))
    order=np.argsort(-sc)
    for i in order:
        s=sc[i]; sx,sy=(s, s*rho) if ori[i]==0 else (s*rho, s)
        ax.add_patch(Ellipse((C[i][0], C[i][1]), 2*sx, 2*sy, angle=0, facecolor=col_i(i),
                             edgecolor=EDGE, lw=min(0.6, 30*s), alpha=0.97))
        if i < a.label_n:                                          # label the first N primes
            fs=min(a.font_max, a.font_base + a.font_scale*s)
            ax.text(C[i][0], C[i][1], str(primes[i]), ha='center', va='center',
                    color='white', fontweight='bold', fontsize=fs, zorder=50)
    ax.add_patch(Ellipse((0,0), 2.0, 2*rho, angle=0, fill=False, ec=EDGE, lw=1.6))
    ax.set_xlim(-1.05,1.05); ax.set_ylim(-rho*1.05, rho*1.05); ax.set_aspect('equal'); ax.axis('off')
    plt.tight_layout(pad=0.3)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight', pad_inches=0.04)
    print(f"wrote {a.out}  ({n} cells, aspect rho={rho:.4f}, labelled first {min(a.label_n,n)} primes)")

if __name__=="__main__":
    main()
