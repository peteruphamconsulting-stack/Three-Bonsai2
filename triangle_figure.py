#!/usr/bin/env python3
"""
triangle_figure.py -- render the certified least-prime-factor TRIANGLE packing (the prime cells
[p], areas d_p) as a flat figure, coloured by scale, with the largest cells labelled by prime.
Companion to triangle_certify.py; reads the same head JSON (pieces = (p, k, a, ax, ay) in the
rational-edge frame, container vertices (-1,0),(1,0),(0,1)) and reconstructs each cell from its
orientation, anchor corner, and anchor.

Usage:
    python3 triangle_figure.py --in triangle_head_n10906.json --out triangle_packing.png --dpi 300
    python3 triangle_figure.py --in triangle_head_n10906.json --label-n 40    # label the first 40 primes
    python3 triangle_figure.py --in triangle_head_n10906.json --label-n 0     # no labels
Requires: numpy, sympy, matplotlib.
"""
import argparse, json
from math import sqrt
import numpy as np, sympy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

RAMP = ["#0d2a33", "#164a52", "#1f6e6a", "#2f8f7d", "#57ab8e", "#8fc9a8"]   # deep -> light
EDGE = "#20302e"

H = sqrt(2.0)/2.0
COS = [1, H, 0, -H, -1, -H, 0, H]
SIN = [0, H, 1, H, 0, -H, -1, -H]
BASE = np.array([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
CONTAINER = np.array([(-1.0, 0.0), (1.0, 0.0), (0.0, 1.0)])

def piece_vertices(k, a, ax, ay, L):
    c, s = COS[k], SIN[k]
    V = BASE @ np.array([[c, -s], [s, c]]).T
    return np.array([ax, ay]) + L*(V - V[a])

def lpf_shares_float(n):
    """d_p for the first n primes (float; drawing only -- rigor lives in the certifier)."""
    hi = int(sympy.prime(n)); a = 1.0; out = []
    for q in sympy.primerange(2, hi+1):
        out.append(a/q); a *= (1.0 - 1.0/q)
    return out[:n]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", required=True)
    ap.add_argument("--out", default="triangle_packing.png"); ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--label-n", type=int, default=40,
                    help="label the cells of the first N primes (0 = no labels)")
    ap.add_argument("--font-base", type=float, default=6.0, help="min label font size (pt)")
    ap.add_argument("--font-scale", type=float, default=42.0, help="label font growth with cell size")
    ap.add_argument("--font-max", type=float, default=40.0, help="max label font size (pt)")
    a = ap.parse_args()
    d = json.load(open(a.infile))
    pieces = d["pieces"]; n = len(pieces)
    primes = [r[0] for r in pieces]
    shares = lpf_shares_float(n)
    sc = np.sqrt(np.array(shares))                     # linear scale sqrt(d_p)
    legs = np.sqrt(2.0*np.array(shares))
    # colour on a log-log scale of the prime (as in the rectangle/ellipse figures)
    ll = np.log(np.log(np.array(primes, float) + 1.0))
    lo, hi = ll.min(), ll.max()
    def col_i(i):
        t = (ll[i]-lo)/(hi-lo+1e-12)
        return RAMP[min(len(RAMP)-1, int(t*len(RAMP)))]
    fig, ax = plt.subplots(figsize=(11, 11*0.56))
    order = np.argsort(-sc)                            # big first, small on top
    for i in order:
        V = piece_vertices(pieces[i][1], pieces[i][2], pieces[i][3], pieces[i][4], legs[i])
        ax.add_patch(Polygon(V, closed=True, facecolor=col_i(i),
                             edgecolor=EDGE, lw=min(0.6, 30*sc[i]), alpha=0.97))
        if i < a.label_n:
            cx, cy = V.mean(axis=0)
            fs = min(a.font_max, a.font_base + a.font_scale*sc[i])
            ax.text(cx, cy, str(primes[i]), ha='center', va='center',
                    color='white', fontweight='bold', fontsize=fs, zorder=50)
    ax.add_patch(Polygon(CONTAINER, closed=True, fill=False, ec=EDGE, lw=1.6, zorder=200))
    ax.set_xlim(-1.03, 1.03); ax.set_ylim(-0.02, 1.03); ax.set_aspect('equal'); ax.axis('off')
    plt.tight_layout(pad=0.3)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight', pad_inches=0.04)
    print(f"wrote {a.out}  ({n} cells, labelled first {min(a.label_n, n)} primes)")

if __name__ == "__main__":
    main()
