#!/usr/bin/env python3
"""
triangle_recursion_figure.py -- render the least-prime-factor TRIANGLE packing and its
self-similar recursion (the "integers in triangular cells").

Each prime cell [p] is drawn from the container head; the big cells 2,3,5,7 (and, where a state
head is supplied, their descendants) are subdivided by mapping the corresponding state-s head
inside them via the EXACT similarity that carries the container onto the cell. Cell [p1..pj] is
labelled by the integer p1*...*pj and coloured by recursion depth.

The similarity is the orientation-preserving map fixed by the vertex correspondence: the
container's right-angle vertex (0,1) goes to the cell's right-angle vertex, and the two leg
endpoints follow in counterclockwise order. Since a state-s cell is a similar copy of the
container and its children carry the state-s areas sigma_s(q)=d_q/R_<s, the same head,
transformed, subdivides every cell of that state (Prop 4.2); cell [2] recurses with the
container (state-2) head itself.

Inputs (JSON heads produced by triangle_build.py):
    --container triangle_head_n10906.json                     (state 2; required)
    --state3 tri_state_3.json --state5 tri_state_5.json --state7 tri_state_7.json
    --max-level 4        recursion depth (1=primes,2=pq,3=..,4=..)
    --recurse-area 0.006 only subdivide cells whose area fraction exceeds this
    --label-area 0.0016  only label cells whose area fraction exceeds this
    --min-draw 3e-7      skip children below this area fraction (sub-pixel)
    --out triangle_recursion.png   --dpi 200

Usage:
    python3 triangle_recursion_figure.py --container triangle_head_n10906.json \
        --state3 tri_state_3.json --state5 tri_state_5.json --state7 tri_state_7.json \
        --out triangle_recursion.png

Requires: numpy, sympy, matplotlib.
"""
import argparse, json
from math import sqrt
import numpy as np, sympy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# depth palette (level 1 primes ... deeper), matching the rectangle/ellipse recursion figures
LEVEL_COLORS = ["#4a9a8f", "#7d6cbc", "#e08a3c", "#5b8bb0", "#c1584e", "#9bbf6a"]
EDGE = "#2c2c2a"

H = sqrt(2.0)/2.0
COS = [1, H, 0, -H, -1, -H, 0, H]
SIN = [0, H, 1, H, 0, -H, -1, -H]
BASE = np.array([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
CONTAINER = np.array([(-1.0, 0.0), (1.0, 0.0), (0.0, 1.0)])
# container vertex correspondence: right angle C3=(0,1), then legs ccw C1=(-1,0), C2=(1,0);
# [C1-C3 | C2-C3] = [[-1,1],[-1,-1]], inverse below (det = 2)
MINV = 0.5*np.array([[-1.0, -1.0], [1.0, -1.0]])
C3 = np.array([0.0, 1.0])

def piece_vertices(k, a, ax, ay, L):
    """(V0,V1,V2): right angle first, counterclockwise (template order is preserved)."""
    c, s = COS[k], SIN[k]
    V = BASE @ np.array([[c, -s], [s, c]]).T
    return np.array([ax, ay]) + L*(V - V[a])

def cell_similarity(V):
    """The orientation-preserving similarity g(z)=Az+t with g(C3)=V0, g(C1)=V1, g(C2)=V2."""
    A = np.column_stack([V[1]-V[0], V[2]-V[0]]) @ MINV
    t = V[0] - A @ C3
    return A, t

def state_shares_float(s, n):
    """sigma_s(q)=d_q/R_<s for the first n primes q>=s (float; drawing only)."""
    base = int(sympy.primepi(s-1)); hi = int(sympy.prime(base+n))
    a = 1.0; out = []
    for q in sympy.primerange(s, hi+1):
        out.append((int(q), a/q)); a *= (1.0 - 1.0/q)
    return out[:n]

def load_head(path):
    d = json.load(open(path))
    s = int(d.get('state', 2)); pieces = d['pieces']; n = len(pieces)
    sh = state_shares_float(s, n)
    cells = []
    for i, (p, k, a, ax, ay) in enumerate(pieces):
        q, share = sh[i]
        assert q == int(p), f"{path}: piece {i} prime {p} != expected {q} (state {s})"
        cells.append(dict(prime=int(p), k=int(k), a=int(a), ax=float(ax), ay=float(ay), share=share))
    return dict(state=s, cells=cells)

def draw(ax, heads, max_level, recurse_area, label_area, min_draw):
    def rec(A, t, S, state, prefix, level):
        for c in heads[state]['cells']:
            rel = (S*S)*c['share']                       # area fraction relative to root
            if rel < min_draw: continue
            L = sqrt(2.0*c['share'])
            V = piece_vertices(c['k'], c['a'], c['ax'], c['ay'], L)   # in the state frame
            W = V @ A.T + t                              # on the root canvas
            integer = prefix*c['prime']
            col = LEVEL_COLORS[min(level-1, len(LEVEL_COLORS)-1)]
            recursed = (level < max_level and c['prime'] in heads and rel > recurse_area)
            ax.add_patch(Polygon(W, closed=True, facecolor=col, edgecolor=EDGE,
                                 lw=min(0.7, 40*S*sqrt(c['share'])),
                                 alpha=(0.16 if recursed else 0.96), zorder=level))
            if recursed:
                Ac, tc = cell_similarity(V)              # state frame -> cell
                rec(A @ Ac, A @ tc + t, S*sqrt(c['share']), c['prime'], integer, level+1)
            elif rel > label_area:
                cx, cy = W.mean(axis=0)
                ax.text(cx, cy, str(integer), ha='center', va='center', zorder=100,
                        fontsize=min(15, 2.6 + 64*S*sqrt(c['share'])), color='white', weight='bold')
    rec(np.eye(2), np.zeros(2), 1.0, 2, 1, 1)
    ax.add_patch(Polygon(CONTAINER, closed=True, fill=False, ec=EDGE, lw=1.6, zorder=200))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--container", required=True)
    ap.add_argument("--state3"); ap.add_argument("--state5"); ap.add_argument("--state7")
    ap.add_argument("--max-level", type=int, default=4)
    ap.add_argument("--recurse-area", type=float, default=0.006)
    ap.add_argument("--label-area", type=float, default=0.0016)
    ap.add_argument("--min-draw", type=float, default=3e-7)
    ap.add_argument("--out", default="triangle_recursion.png"); ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    heads = {2: load_head(a.container)}
    for s, f in [(3, a.state3), (5, a.state5), (7, a.state7)]:
        if f: heads[s] = load_head(f)
    fig, ax = plt.subplots(figsize=(11, 11*0.60))
    draw(ax, heads, a.max_level, a.recurse_area, a.label_area, a.min_draw)
    ax.set_xlim(-1.03, 1.03); ax.set_ylim(-0.02, 1.06); ax.set_aspect('equal'); ax.axis('off')
    lv = min(a.max_level, 4)
    handles = [plt.Line2D([0], [0], marker='^', ls='', mfc=LEVEL_COLORS[i], mec=EDGE, ms=11) for i in range(lv)]
    labels = ["primes $p$", "$p\\cdot q$", "3 factors", "4 factors"][:lv]
    ax.legend(handles, labels, loc='upper right', frameon=False, fontsize=11, title="recursion depth")
    ax.set_title("The integers in the triangle: recursive least-prime-factor packing", fontsize=13)
    plt.tight_layout(pad=0.4)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight', pad_inches=0.05)
    print(f"wrote {a.out}  (states with heads: {sorted(heads)})")

if __name__ == "__main__":
    main()
