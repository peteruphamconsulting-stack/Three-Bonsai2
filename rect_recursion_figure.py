#!/usr/bin/env python3
"""
rect_recursion_figure.py -- the natural numbers in the sqrt(2) rectangle: the recursive
least-prime-factor guillotine dissection ("integers as a point set").

Self-contained (needs only a prime sieve): the guillotine packing is exact and deterministic,
so no head files are required. The container sqrt(2)-rectangle is guillotine-packed by the prime
cells [p] (areas d_p); each cell of state s (largest prime factor) is a sqrt(2)-rectangle and is
packed by the SAME guillotine on its children's relative areas sigma_s(q)=d_q/R_<s -- exact
self-similarity (Prop 4.2). Cells are labelled by their integer p1*...*pj and coloured by depth.

    area[n] = R_<P(n) / n  (P(n) = largest prime factor): powers of two are the largest composite
    cells, primes the most discounted -- ordered by largest prime factor, not by n.

Usage:
    python3 rect_recursion_figure.py --out rect_recursion.png --dpi 250 \
        --max-level 4 --recurse-area 0.004 --min-rel 0.004 --label-area 0.0016
Requires: numpy, sympy, matplotlib.
"""
import argparse
from math import sqrt
import numpy as np, sympy
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SQ2 = sqrt(2.0)
LEVEL_COLORS = ["#4a9a8f", "#7d6cbc", "#e08a3c", "#5b8bb0", "#c1584e", "#9bbf6a"]
RESID = "#eceae3"
EPS = 1e-13

def state_children(s, min_rel):
    """(q, sigma_s(q)) for primes q>=s with sigma_s(q)>=min_rel; sigma_s decreasing in q.
    sigma_s(q) = (1/q) prod_{s<=r<q}(1-1/r) = d_q / R_<s."""
    out=[]; a=1.0
    for q in sympy.primerange(s, 10**7):
        rel=a/q
        if rel < min_rel and q > s: break
        out.append((int(q), rel)); a *= (1.0 - 1.0/q)
    return out

def guillotine_pack(rect, areas, r=SQ2):
    """Pack sqrt(2)-rectangles of the given areas (absolute) into rect=(x,y,W,H) by the shelf
    guillotine of fig_rect. Returns (placed=[(x,y,pw,ph,idx)], free=[(x,y,w,h)])."""
    x0,y0,W,H = rect
    fr=[(x0,y0,W,H)]; placed=[]
    for k,d in enumerate(areas):
        a,b = sqrt(d*r), sqrt(d/r)
        best=None
        for i,(x,y,w,h) in enumerate(fr):
            if b<=min(w,h)+1e-12 and a<=max(w,h)+1e-12:
                if best is None or min(w,h)>min(fr[best][2],fr[best][3]): best=i
        if best is None: continue                        # this piece doesn't fit; skip (residual)
        x,y,w,h = fr.pop(best)
        pw,ph = (a,b) if (a<=w+1e-12 and b<=h+1e-12) else (b,a)
        placed.append((x,y,pw,ph,k))
        c1=max(min(w-pw,h),min(pw,h-ph)); c2=max(min(w-pw,ph),min(w,h-ph))
        if c1>=c2:
            if w-pw>EPS: fr.append((x+pw,y,w-pw,h))
            if h-ph>EPS: fr.append((x,y+ph,pw,h-ph))
        else:
            if w-pw>EPS: fr.append((x+pw,y,w-pw,ph))
            if h-ph>EPS: fr.append((x,y+ph,w,h-ph))
    return placed, fr

def main():
    ap=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="rect_recursion.png"); ap.add_argument("--dpi", type=int, default=250)
    ap.add_argument("--max-level", type=int, default=4)
    ap.add_argument("--recurse-states", default="2,3,5,7")
    ap.add_argument("--recurse-area", type=float, default=0.004)
    ap.add_argument("--min-rel", type=float, default=0.004)
    ap.add_argument("--label-area", type=float, default=0.0016)
    ap.add_argument("--font-base", type=float, default=4.0, help="min label font size (pt)")
    ap.add_argument("--font-scale", type=float, default=26.0, help="label font growth with cell size")
    ap.add_argument("--font-max", type=float, default=22.0, help="max label font size (pt)")
    a=ap.parse_args()
    rstates=set(int(s) for s in a.recurse_states.split(","))
    W,H = sqrt(SQ2), 1/sqrt(SQ2)                          # container sqrt(2)-rectangle, area 1
    Aroot = W*H
    fig,ax=plt.subplots(figsize=(10, 10*H/W))

    def rec(rect, state, prefix, level):
        kids = state_children(state, a.min_rel)
        cellA = rect[2]*rect[3]
        placed, free = guillotine_pack(rect, [rel*cellA for _,rel in kids])
        for (x,y,pw,ph,k) in placed:
            q = kids[k][0]; integer = prefix*q
            rel = (pw*ph)/Aroot
            col = LEVEL_COLORS[min(level-1, len(LEVEL_COLORS)-1)]
            recursed = (level < a.max_level and q in rstates and rel > a.recurse_area)
            ax.add_patch(Rectangle((x,y), pw, ph, facecolor=col, edgecolor="white",
                                   lw=0.4, alpha=(0.16 if recursed else 0.96), zorder=level))
            if recursed:
                rec((x,y,pw,ph), q, integer, level+1)
            elif rel > a.label_area:
                fs=min(a.font_max, a.font_base + a.font_scale*sqrt(rel))
                ax.text(x+pw/2, y+ph/2, str(integer), ha="center", va="center", zorder=100,
                        color="white", fontweight="bold", fontsize=fs)
        for (x,y,w,h) in free:                            # residual (smaller children, not drawn)
            ax.add_patch(Rectangle((x,y), w, h, facecolor=RESID, edgecolor="none", zorder=0))

    rec((0.0,0.0,W,H), 2, 1, 1)
    ax.add_patch(Rectangle((0,0), W, H, fill=False, edgecolor="k", lw=1.5, zorder=200))
    ax.set_xlim(-0.01, W+0.01); ax.set_ylim(-0.01, H+0.01); ax.set_aspect("equal"); ax.axis("off")
    handles=[plt.Line2D([0],[0],marker='s',ls='',mfc=LEVEL_COLORS[i],mec='white',ms=12)
             for i in range(min(a.max_level,4))]
    labels=["primes $p$","$p\\cdot q$","3 factors","4 factors"][:min(a.max_level,4)]
    leg=ax.legend(handles,labels,loc='center left',bbox_to_anchor=(1.01, 0.5),
                  frameon=False,fontsize=11,title="recursion depth")
    ax.set_title(r"The integers in the $\sqrt{2}$ rectangle: recursive least-prime-factor guillotine",
                 fontsize=13)
    plt.tight_layout(pad=0.4)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches="tight", pad_inches=0.05, bbox_extra_artists=(leg,))
    print(f"wrote {a.out}")

if __name__=="__main__":
    main()
