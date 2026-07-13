#!/usr/bin/env python3
"""
ellipse_recursion_figure.py -- render the least-prime-factor ELLIPSE packing and its
self-similar recursion (the "integers as a point set").

Each prime cell [p] is drawn from the certified container head; the big cells 2,3,5,7 (and,
where a state head is supplied, their descendants) are subdivided by mapping the corresponding
state-s head inside them via the exact similarity that carries the unit container onto the cell.
Cell [p1..pj] is labelled by the integer p1*...*pj and coloured by recursion depth.

The self-similarity is exact: a cell of state s is a scaled/rotated copy of the aspect-rho
container, and its children carry the state-s areas sigma_s(q)=d_q/R_<s; so the same head,
transformed, tiles every cell of that state.

Inputs (JSON heads produced by ellipse_build.py):
    --container ellipse_head_n1650.json          (state 2; required)
    --state3 head_s3.json  --state5 head_s5.json  --state7 head_s7.json   (optional; deeper layers)
    --max-level 4        recursion depth (1=primes,2=pq,3=..,4=..)
    --recurse-area 0.008 only subdivide cells whose area fraction exceeds this
    --label-area 0.0016  only label cells whose area fraction exceeds this
    --out ellipse_recursion.png   --dpi 200

Usage:
    python3 ellipse_recursion_figure.py --container ellipse_head_n1650.json \
        --state3 head_s3.json --state5 head_s5.json --state7 head_s7.json --out ellipse_recursion.png

Requires: numpy, sympy, matplotlib.
"""
import argparse, json
from math import sqrt, pi
import numpy as np, sympy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# depth palette (level 1 primes ... deeper), tuned to read like the LPF-tree figures
LEVEL_COLORS = ["#4a9a8f", "#7d6cbc", "#e08a3c", "#5b8bb0", "#c1584e", "#9bbf6a"]
EDGE = "#2c2c2a"

def first_primes_from(s, n):
    base=int(sympy.primepi(s-1)); hi=int(sympy.prime(base+n))
    return [int(p) for p in sympy.primerange(s, hi+1)][:n]

def load_head(path):
    d=json.load(open(path)); V=d.get('verification', d)
    rho=float(V['aspect_ratio_b_over_a']); s=int(V.get('state', d.get('state',2)))
    ori=[1 if o==90 else 0 for o in V['orientation_pattern_degrees']]
    shares=V['area_shares_float']; n=len(shares)
    primes=first_primes_from(s, n)
    cells=[]
    for i in range(n):
        cx,cu = V['centers_normalized'][i] if 'centers_normalized' in V else (V['centers'][i][0], V['centers'][i][1]/rho)
        cells.append(dict(x=cx, u=cu, ori=ori[i], prime=primes[i], share=shares[i]))
    return dict(rho=rho, state=s, cells=cells)

def draw(ax, heads, root_state, max_level, recurse_area, label_area):
    rho = heads[root_state]['rho']
    # Work in ACTUAL Euclidean coords. A cell of state s is the unit container scaled by S,
    # rotated by R in {0,90} deg, centred at (CX,CY). Mapping the state-s head inside it is a
    # genuine Euclidean similarity: rotate the child's actual centre, swap its semi-axes if R=90.
    def rec(CX, CY, S, R, state, prefix, level):
        for c in heads[state]['cells']:
            X0, Y0 = c['x'], rho*c['u']                 # child actual centre in the unit-container frame
            s0 = sqrt(c['share'])
            sx0, sy0 = (s0, s0*rho) if c['ori']==0 else (s0*rho, s0)   # child actual semi-axes
            if R == 0:
                X, Y = CX + S*X0, CY + S*Y0
                SX, SY, nR = S*sx0, S*sy0, c['ori']
            else:                                        # 90 deg CCW: (X0,Y0)->(-Y0,X0); axes swap
                X, Y = CX - S*Y0, CY + S*X0
                SX, SY, nR = S*sy0, S*sx0, 1 - c['ori']
            integer = prefix * c['prime']
            rel = (S*s0)**2                              # area fraction relative to root container
            col = LEVEL_COLORS[min(level-1, len(LEVEL_COLORS)-1)]
            recursed = (level < max_level and c['prime'] in heads and rel > recurse_area)
            ax.add_patch(Ellipse((X, Y), 2*SX, 2*SY, angle=0, facecolor=col, edgecolor=EDGE,
                                 lw=min(0.7, 40*S*s0), alpha=(0.16 if recursed else 0.96), zorder=level))
            if recursed:
                rec(X, Y, S*s0, nR, c['prime'], integer, level+1)
            elif rel > label_area:
                ax.text(X, Y, str(integer), ha='center', va='center', zorder=100,
                        fontsize=min(15, 2.6+90*S*s0), color='white', weight='bold')
    rec(0.0, 0.0, 1.0, 0, root_state, 1, 1)
    ax.add_patch(Ellipse((0,0), 2*1.0, 2*rho, angle=0, fill=False, ec=EDGE, lw=1.6, zorder=200))

def main():
    ap=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--container", required=True)
    ap.add_argument("--state3"); ap.add_argument("--state5"); ap.add_argument("--state7")
    ap.add_argument("--max-level", type=int, default=4)
    ap.add_argument("--recurse-area", type=float, default=0.008)
    ap.add_argument("--label-area", type=float, default=0.0016)
    ap.add_argument("--out", default="ellipse_recursion.png"); ap.add_argument("--dpi", type=int, default=200)
    a=ap.parse_args()
    heads={2: load_head(a.container)}
    for s,f in [(3,a.state3),(5,a.state5),(7,a.state7)]:
        if f: heads[s]=load_head(f)
    rho=heads[2]['rho']
    fig,ax=plt.subplots(figsize=(11, 11*rho*1.08))
    draw(ax, heads, 2, a.max_level, a.recurse_area, a.label_area)
    ax.set_xlim(-1.05,1.05); ax.set_ylim(-rho*1.05, rho*1.05); ax.set_aspect('equal'); ax.axis('off')
    lv=[f"level {i+1}" for i in range(min(a.max_level,4))]
    handles=[plt.Line2D([0],[0],marker='o',ls='',mfc=LEVEL_COLORS[i],mec=EDGE,ms=11) for i in range(len(lv))]
    labels=["primes $p$","$p\\cdot q$","3 factors","4 factors"][:len(lv)]
    ax.legend(handles,labels,loc='upper right',frameon=False,fontsize=11,title="recursion depth")
    ax.set_title("The integers in the ellipse: recursive least-prime-factor packing", fontsize=13)
    plt.tight_layout(pad=0.4)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight', pad_inches=0.05)
    print(f"wrote {a.out}  (states with heads: {sorted(heads)})")

if __name__=="__main__":
    main()
