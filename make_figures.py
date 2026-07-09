#!/usr/bin/env python3
"""
make_figures.py  --  regenerate every figure in
"A Geometric Sieve of Eratosthenes" (P. Upham).

Figures produced (into the current directory):
    staircase_angle.png   (Fig 1, Thm 3.1 : limit angle of the staircase model)
    integers.png          (Fig 2       : three levels of the integer dissection)
    packing.png           (Fig 3       : 355-piece corner-anchored triangle greedy)
    chain.png             (Fig 4       : 39-piece chain packing + centroid path)
    certificate.png       (Fig 5       : free-region decomposition + supply inequality)
    rect_packing.png      (Fig 6       : exact guillotine packing of the sqrt2 rectangle)
    vertex_stats.png      (Fig 7       : vertex statistics of the 355-piece packing)

Data files expected in the same directory (produced by the companion enumeration code):
    pieces_greedy.json      list of {poly:[3 pts], p, d}   -- 355-piece triangle greedy
    pieces_chain_dfs.json   list of {poly:[3 pts], p, d}   -- 39-piece chain
    integers.json           list of {poly, label, level, kind}
    head_certificate.json   dict with regions/widths/region_areas (p_last=2389)
The staircase and rectangle figures are self-contained (need only a prime sieve).

Usage:  python3 make_figures.py [figure_name ...]
        (no args -> all figures)
"""
import sys, json
from math import sqrt, log, atan2, pi, exp, degrees
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon as MplPoly, Rectangle

SQ2 = sqrt(2.0)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)

def lpf_densities(P):
    """d_p = (1/p) prod_{q<p}(1-1/q)  for the primes in array P."""
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0 / P.astype(float)))[:-1]])
    return np.exp(logfac) / P, np.exp(logfac)      # d_p , R_{<p}

def load(name):
    return json.load(open(name))

def free_region(polys, container):
    """container minus union of triangle polygons, via shapely (float)."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    T = Polygon(container)
    return T.difference(unary_union([Polygon(p) for p in polys]))


# --------------------------------------------------------------------------
# Fig 1 : staircase limit angle (self-contained)
# --------------------------------------------------------------------------
def fig_staircase(path="staircase_angle.png", N=48):
    P = sieve(10**4)[:N]
    # remaining rectangle anchored at top-right corner M=(1,1); alternate cuts.
    # odd step j: vertical cut -> w *= (1-1/p_j);  even step: horizontal -> h *= (1-1/p_j)
    w = h = 1.0
    xs, ys, ratios = [1.0 - w], [1.0 - h], []
    for j, p in enumerate(P, start=1):
        if j % 2 == 1: w *= (1.0 - 1.0 / p)
        else:          h *= (1.0 - 1.0 / p)
        xs.append(1.0 - w); ys.append(1.0 - h); ratios.append(h / w)
    # the constant L (alternating prime sum) converges slowly: use many primes for its value
    Pbig = sieve(2 * 10**7)
    L = float(np.sum(((-1.0)**np.arange(1, len(Pbig) + 1)) * np.log(1.0 - 1.0 / Pbig)))
    kappa = exp(L); theta = degrees(atan2(kappa, 1.0))

    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    # left: staircase vertices + limiting ray
    ax[0].step(xs, ys, where="post", color="0.6", lw=0.8)
    ax[0].scatter(xs, ys, s=9, c=range(len(xs)), cmap="viridis")
    t = np.linspace(0, 1, 2)
    ax[0].plot(1 - t * (1 / sqrt(1 + kappa**2)), 1 - t * (kappa / sqrt(1 + kappa**2)),
               "r--", lw=1.2, label=fr"limiting ray $\theta^*={theta:.2f}^\circ$")
    ax[0].scatter([1], [1], marker="*", s=120, color="k", zorder=5, label="limit point $M$")
    ax[0].set_aspect("equal"); ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
    ax[0].legend(loc="lower left", fontsize=9); ax[0].set_title("staircase vertices")
    # right: aspect ratio -> kappa
    ax[1].axhline(kappa, color="r", lw=1, label=fr"$\kappa=e^{{L}}={kappa:.5f}$")
    ax[1].plot(range(1, len(ratios) + 1), ratios, ".-", color="teal", ms=4)
    ax[1].set_xscale("log"); ax[1].set_xlabel("$n$"); ax[1].set_ylabel("$h_n/w_n$")
    ax[1].legend(fontsize=9); ax[1].set_title("convergence of the aspect ratio")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: L={L:.10f}, kappa={kappa:.7f}, theta*={theta:.3f} deg")


# --------------------------------------------------------------------------
# Fig 2 : integer dissection (three levels)
# --------------------------------------------------------------------------
def fig_integers(path="integers.png"):
    cells = load("integers.json")
    lvl_color = {1: "#41ab5d", 2: "#807dba", 3: "#fd8d3c", 4: "#fdd0a2"}   # green / violet / gold
    def area(poly):
        p = np.array(poly); return 0.5*abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]))
    def leg(poly):                                   # short side ~ cell scale
        p = np.array(poly); e = [np.linalg.norm(p[(i+1)%3]-p[i]) for i in range(3)]
        return min(e)
    def ra_vertex(poly):                             # right-angle corner (bracket anchor)
        p = np.array(poly); best=0; bv=9
        for i in range(3):
            u=p[(i+1)%3]-p[i]; v=p[(i-1)%3]-p[i]
            c=abs(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)+1e-12))
            if c<bv: bv=c; best=i
        return p[best]
    fig, ax = plt.subplots(figsize=(9, 9))
    for c in sorted(cells, key=lambda c: -area(c["poly"])):     # big first, small on top
        poly = np.array(c["poly"]); fc = lvl_color.get(c["level"], "0.8")
        ax.add_patch(MplPoly(poly, closed=True, facecolor=fc, edgecolor="white", lw=0.5,
                             alpha=0.42 if c["kind"] == "shell" else 0.93))
    boxes = []  # (x,y,hw,hh) of drawn labels, to keep brackets clear of numbers
    s = SQ2 / 8.0
    def collides(x, y, hw, hh):
        return any(abs(x-bx) < hw+bw and abs(y-by) < hh+bh for bx,by,bw,bh in boxes)
    # leaves first (numbers), gated by cell size
    for c in cells:
        if c["kind"] != "leaf": continue
        poly = np.array(c["poly"]); L = leg(poly)
        if L <= 0.055: continue
        fs = float(np.clip(L*46, 6, 15)); cen = poly.mean(0)
        hw = 0.5*len(str(c["label"]))*fs*0.6/72*s; hh = 0.5*fs/72*s
        ax.text(cen[0], cen[1], str(c["label"]), ha="center", va="center",
                fontsize=fs, fontweight="bold", color="white")
        boxes.append((cen[0], cen[1], hw, hh))
    # shell brackets, placed only where clear of numbers
    for c in cells:
        if c["kind"] != "shell": continue
        poly = np.array(c["poly"]); L = leg(poly)
        if L <= 0.16: continue
        v = ra_vertex(poly); cen = poly.mean(0); d = cen - v; d = d/(np.linalg.norm(d)+1e-9)
        placed = False
        for f in (0.10, 0.16, 0.24, 0.34, 0.46):     # walk in from the tip toward centroid
            x, y = v[0]+f*(cen[0]-v[0])/1, v[1]+f*(cen[1]-v[1])/1
            if not collides(x, y, 0.028, 0.018):
                ax.text(x, y, f"[{c['label']}]", fontsize=7.5, color="0.10",
                        ha="center", va="center", fontweight="bold")
                boxes.append((x, y, 0.028, 0.018)); placed = True; break
        # if nowhere clear, leave the shell unlabeled (its colour/transparency already marks it)
    ax.add_patch(MplPoly([(0, 0), (SQ2, 0), (0, SQ2)], closed=True, fill=False, edgecolor="k", lw=1.3))
    ax.set_xlim(-0.02, SQ2 + 0.02); ax.set_ylim(-0.02, SQ2 + 0.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("The natural numbers in the triangle: recursive least-prime-factor sieve\\n"
                 "level 1 (green) $n=p$ · level 2 (violet) $n=pq$ · level 3 (gold) $n=4q$", fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {len(cells)} cells, "
          f"{sum(c['kind']=='shell' for c in cells)} shells, {sum(c['kind']=='leaf' for c in cells)} leaves")


# --------------------------------------------------------------------------
# Fig 3 : triangle greedy packing
# --------------------------------------------------------------------------
def fig_packing(path="packing.png"):
    pieces = load("pieces_greedy.json")
    fig, ax = plt.subplots(figsize=(8, 8))
    lp = max(q["p"] for q in pieces)
    for q in pieces:
        poly = np.array(q["poly"])
        col = cm.viridis((log(log(q["p"] + 1)) - log(log(2))) /
                         (log(log(lp)) - log(log(2))))
        ax.add_patch(MplPoly(poly, closed=True, facecolor=col,
                             edgecolor="white", lw=0.3))
        if q["p"] <= 23:
            ax.text(*poly.mean(0), str(q["p"]), ha="center", va="center",
                    color="white", fontsize=12, fontweight="bold")
    free = free_region([q["poly"] for q in pieces], [(0, 0), (SQ2, 0), (0, SQ2)])
    for g in (free.geoms if free.geom_type != "Polygon" else [free]):
        if g.area > 1e-13:
            ax.add_patch(MplPoly(np.array(g.exterior.coords), closed=True,
                                 facecolor="#f4a582", edgecolor="none"))
    ax.set_xlim(-0.02, SQ2 + 0.02); ax.set_ylim(-0.02, SQ2 + 0.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Corner-anchored greedy: {len(pieces)} pieces, primes 2–{lp}, "
                 f"free region (orange) one component")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {len(pieces)} pieces to p={lp}")


# --------------------------------------------------------------------------
# Fig 4 : chain packing + centroid path
# --------------------------------------------------------------------------
def fig_chain(path="chain.png"):
    pieces = load("pieces_chain_dfs.json")
    fig, ax = plt.subplots(figsize=(8, 8))
    cents = []
    lp = max(q["p"] for q in pieces)
    for q in pieces:
        poly = np.array(q["poly"]); cents.append(poly.mean(0))
        col = cm.plasma((log(q["p"]) - log(2)) / (log(lp) - log(2)))
        ax.add_patch(MplPoly(poly, closed=True, facecolor=col,
                             edgecolor="white", lw=0.4))
    cents = np.array(cents)
    ax.plot(cents[:, 0], cents[:, 1], "-", color="k", lw=1.0, alpha=0.7)
    ax.scatter(cents[:, 0], cents[:, 1], s=10, color="k", zorder=5)
    ax.add_patch(MplPoly([(0, 0), (SQ2, 0), (0, SQ2)], closed=True,
                         fill=False, edgecolor="k", lw=1.2))
    ax.set_xlim(-0.02, SQ2 + 0.02); ax.set_ylim(-0.02, SQ2 + 0.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Chain packing: {len(pieces)} pieces, primes 2–{lp}; centroid path drawn")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {len(pieces)} pieces to p={lp}")


# --------------------------------------------------------------------------
# Fig 5 : free-region decomposition + supply inequality
# --------------------------------------------------------------------------
def fig_certificate(path="certificate.png"):
    H = load("head_certificate.json")
    regions = H["regions"]; widths = np.array(H["widths"])
    areas = np.array(H["region_areas"]); p_last = H["p_last"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    # left: regions coloured by width (log scale)
    wln = np.log10(widths + 1e-12)
    norm = (wln - wln.min()) / (wln.max() - wln.min() + 1e-12)
    for poly, u in zip(regions, norm):
        ax[0].add_patch(MplPoly(np.array(poly), closed=True,
                                facecolor=cm.viridis(u), edgecolor="white", lw=0.15))
    ax[0].add_patch(MplPoly([(0, 0), (SQ2, 0), (0, SQ2)], closed=True,
                            fill=False, edgecolor="k", lw=1.0))
    ax[0].set_xlim(-0.02, SQ2 + 0.02); ax[0].set_ylim(-0.02, SQ2 + 0.02)
    ax[0].set_aspect("equal"); ax[0].axis("off")
    ax[0].set_title(f"free region at $p_0={p_last}$: {len(regions)} convex cells,\ncoloured by width")
    # right: supply inequality sigma_<=(w) <= R(p(w))
    P = sieve(3 * 10**6); dP, Rlt = lpf_densities(P)
    order = np.argsort(widths)
    w_sorted = widths[order]; cum_area = np.cumsum(areas[order])
    # future supply: pieces small enough to enter width w have a_p<=w i.e. sqrt(d_p)<=w
    aP = np.sqrt(dP)
    R_supply = np.interp(w_sorted, aP[::-1], Rlt[::-1])   # R at the prime with a_p=w
    ax[1].loglog(w_sorted, cum_area, color="crimson", lw=1.6,
                 label=r"$\sigma_{\leq}(w)$  (free area in cells of width $\leq w$)")
    ax[1].loglog(w_sorted, R_supply, color="navy", lw=1.2, ls="--",
                 label=r"$R_{<p(w)}$  (future supply)")
    ax[1].set_xlabel("width $w$"); ax[1].set_ylabel("area")
    ax[1].legend(fontsize=9); ax[1].set_title("supply inequality")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {len(regions)} regions at p_0={p_last}")


# --------------------------------------------------------------------------
# Fig 6 : exact guillotine packing of the sqrt2 rectangle (self-contained)
# --------------------------------------------------------------------------
def fig_rect(path="rect_packing.png", N=430, r=SQ2):
    import heapq
    P = sieve(2 * 10**5); dP, _ = lpf_densities(P)
    fr = [(0.0, 0.0, sqrt(r), 1 / sqrt(r))]      # (x, y, w, h)
    placed = []
    for k in range(N):
        d = dP[k]; a, b = sqrt(d * r), sqrt(d / r)
        best = None
        for i, (x, y, w, h) in enumerate(fr):
            if b <= min(w, h) + 1e-15 and a <= max(w, h) + 1e-15:
                if best is None or min(w, h) > min(fr[best][2], fr[best][3]):
                    best = i
        x, y, w, h = fr.pop(best)
        pw, ph = (a, b) if (a <= w + 1e-15 and b <= h + 1e-15) else (b, a)
        placed.append((x, y, pw, ph, int(P[k])))
        c1 = max(min(w - pw, h), min(pw, h - ph))
        c2 = max(min(w - pw, ph), min(w, h - ph))
        if c1 >= c2:
            if w - pw > 1e-13: fr.append((x + pw, y, w - pw, h))
            if h - ph > 1e-13: fr.append((x, y + ph, pw, h - ph))
        else:
            if w - pw > 1e-13: fr.append((x + pw, y, w - pw, ph))
            if h - ph > 1e-13: fr.append((x, y + ph, w, h - ph))
    fig, ax = plt.subplots(figsize=(9, 9 / r))
    lp = placed[-1][4]
    for (x, y, pw, ph, p) in placed:
        col = cm.viridis((log(log(p + 1)) - log(log(2))) / (log(log(lp)) - log(log(2))))
        ax.add_patch(Rectangle((x, y), pw, ph, facecolor=col, edgecolor="white", lw=0.3))
        if p <= 47:
            ax.text(x + pw / 2, y + ph / 2, str(p), ha="center", va="center",
                    color="white", fontsize=13, fontweight="bold")
    for (x, y, w, h) in fr:
        ax.add_patch(Rectangle((x, y), w, h, facecolor="#f4a582", edgecolor="none"))
    ax.set_xlim(0, sqrt(r)); ax.set_ylim(0, 1 / sqrt(r)); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(fr"Guillotine packing of similar rectangles, areas $d_p$, $r=\sqrt{{2}}$"
                 "\n"
                 fr"{N} pieces, primes 2–{lp} · free region (orange) exact rectangles", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {N} pieces, {len(fr)} free rectangles")


# --------------------------------------------------------------------------
# Fig 7 : vertex statistics (4 panels)
# --------------------------------------------------------------------------
def fig_vertex_stats(path="vertex_stats.png"):
    pieces = load("pieces_greedy.json")
    cents = np.array([np.array(q["poly"]).mean(0) for q in pieces])
    ps = np.array([q["p"] for q in pieces])
    M = cents[-1]                                   # accumulation locus proxy
    d = cents - M; dist = np.hypot(d[:, 0], d[:, 1]); ang = np.arctan2(d[:, 1], d[:, 0])
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].scatter(ps, dist, s=7, c="teal"); ax[0].set_xscale("log")
    ax[0].set_xlabel("prime $p$"); ax[0].set_ylabel("distance to $M$")
    ax[0].set_title("piece distance from accumulation locus")
    ax[1].scatter(ps, ang, s=7, c="indigo"); ax[1].set_xscale("log")
    ax[1].set_xlabel("prime $p$"); ax[1].set_ylabel("angle about $M$")
    ax[1].set_title("angular coordinate")
    ax[2].hist(ang, bins=40, density=True, color="mediumseagreen")
    ax[2].axhline(1 / (2 * pi), color="r", lw=1, label="uniform")
    ax[2].legend(fontsize=9); ax[2].set_xlabel("angle about $M$")
    ax[2].set_title(f"angle distribution, {len(ps)} vertices")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"  {path}: {len(pieces)} vertices (greedy packing)")


FIGS = {
    "staircase_angle": fig_staircase, "integers": fig_integers,
    "packing": fig_packing, "chain": fig_chain, "certificate": fig_certificate,
    "rect_packing": fig_rect, "vertex_stats": fig_vertex_stats,
}

if __name__ == "__main__":
    want = sys.argv[1:] or list(FIGS)
    for name in want:
        fn = FIGS.get(name.replace(".png", ""))
        if fn is None:
            print(f"unknown figure: {name}; choices: {list(FIGS)}"); continue
        print(f"generating {name} ...")
        try:
            fn()
        except FileNotFoundError as e:
            print(f"  SKIPPED ({name}): missing data file {e.filename}")
    print("done.")
