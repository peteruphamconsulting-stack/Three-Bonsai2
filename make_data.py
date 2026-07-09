#!/usr/bin/env python3
"""
make_data.py  --  enumeration companion for the figures of
"A Geometric Sieve of Eratosthenes" (P. Upham).

Produces the JSON data files consumed by make_figures.py:

    pieces_greedy.json     list of {poly:[3 pts], p, d}   -- corner-anchored triangle greedy
    pieces_chain_dfs.json  list of {poly:[3 pts], p, d}   -- adjacency chain (backtracking DFS)
    integers.json          list of {poly, label, level, kind}  -- recursive integer dissection
    head_certificate.json  {regions, widths, region_areas, p_last}  -- free-region decomposition

The two self-contained figures (staircase_angle.png, rect_packing.png) need no data file.
Only integers.json, pieces_greedy.json and pieces_chain_dfs.json are required for the four
figures the current paper uses; head_certificate.json is for the optional certificate figure.

Geometry.  The container is the isosceles right triangle T = (0,0),(sqrt2,0),(0,sqrt2) of area 1.
The piece for prime p is a similar right isosceles triangle of area d_p, i.e. leg L_p = sqrt(2 d_p)
and altitude a_p = sqrt(d_p). Placements are corner-anchored in the 8 axis-aligned 45-degree
orientations; the free region is maintained exactly as a shapely polygon and each piece is chosen
to nestle deepest (longest shared boundary) among all valid anchored placements.

Requires: numpy, shapely (`pip install shapely`).

Usage:
    python3 make_data.py                       # all four files, default sizes
    python3 make_data.py greedy --nmax 355     # just the greedy packing (primes -> ~2389)
    python3 make_data.py chain  --nmax 60
    python3 make_data.py integers
    python3 make_data.py certificate --nmax 355
Options:
    --nmax N        number of primes to place (greedy/chain/certificate)
    --time-cap S    stop the greedy/chain after S seconds (0 = no cap)
Notes:
    The greedy is O(pieces x free-vertices x orientations) shapely containment tests; expect a few
    minutes for ~355 pieces. Counts (355 pieces, 39-piece chain) depend on the exact placement rule
    and tolerances; if your run differs, update the captions in the .tex to the numbers printed here.
"""
import sys, json, argparse, time
from math import sqrt, atan2, pi
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.prepared import prep
from shapely.ops import unary_union

SQ2 = sqrt(2.0)
TOL = 1e-9
UNIT_T = [(0.0, 0.0), (SQ2, 0.0), (0.0, SQ2)]          # right angle at vertex 0


# ---------------------------------------------------------------- primes / densities
def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)

def lpf_densities(P):
    """d_p = (1/p) prod_{q<p}(1-1/q); also return R_{<p} = prod_{q<p}(1-1/q)."""
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0 / P.astype(float)))[:-1]])
    Rlt = np.exp(logfac)
    return Rlt / P, Rlt


# ---------------------------------------------------------------- piece geometry
def orient_variants(L):
    """The distinct corner-anchored placements of a leg-L right isosceles piece:
    8 rotations (multiples of 45 deg) x 3 anchor corners, deduplicated. Each entry is a
    (3,2) array of vertices with the anchor corner translated to the origin, plus the index
    of the right-angle vertex (needed for the affine map in the integer dissection)."""
    base = np.array([(0.0, 0.0), (L, 0.0), (0.0, L)])   # right angle at index 0
    out = []
    for k in range(8):
        th = k * pi / 4; c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s], [s, c]])
        v = base @ R.T
        for a in range(3):                              # anchor corner a at origin
            shifted = v - v[a]
            ra = 0                                      # right-angle vertex is base index 0
            out.append((shifted, ra))
    uniq = []
    for o, ra in out:
        if not any(np.allclose(o, u, atol=1e-12) for u, _ in uniq):
            uniq.append((o, ra))
    return uniq


def free_vertices(F):
    vs = []
    geoms = F.geoms if isinstance(F, MultiPolygon) else [F]
    for g in geoms:
        vs += list(g.exterior.coords)[:-1]
        for r in g.interiors:
            vs += list(r.coords)[:-1]
    return vs


def _clean(F):
    if F.geom_type == 'GeometryCollection':
        F = unary_union([g for g in F.geoms if 'Polygon' in g.geom_type])
    return F


# ---------------------------------------------------------------- greedy packing
def greedy_pack(container, areas, time_cap=0, verbose=True):
    """Corner-anchored deepest-nestle greedy. Returns list of (poly(3x2 array), ra_index)."""
    F = Polygon(container); placed = []; cache = {}; t0 = time.time()
    for idx, area in enumerate(areas):
        L = sqrt(2 * area); key = round(L, 12)
        os = cache.get(key)
        if os is None: os = cache[key] = orient_variants(L)
        pf = prep(F); best = None; bestsc = -1.0
        for v in free_vertices(F):
            vv = np.array(v)
            for o, ra in os:
                poly = Polygon(o + vv)
                if not poly.is_valid: continue
                if pf.contains(poly.buffer(-TOL)):
                    inter = poly.boundary.intersection(F.boundary)
                    sc = inter.length if not inter.is_empty else 0.0
                    if sc > bestsc + 1e-12:
                        bestsc = sc; best = (o + vv, ra, poly)
        if best is None:
            if verbose: print(f"  greedy jam at piece {idx}")
            break
        pts, ra, poly = best
        placed.append((pts, ra))
        F = _clean(F.difference(poly))
        if verbose and (idx + 1) % 25 == 0:
            print(f"  {idx+1} pieces, covered {1-F.area:.4f}, {time.time()-t0:.0f}s")
        if time_cap and time.time() - t0 > time_cap:
            if verbose: print(f"  time cap after {idx+1} pieces"); break
    return placed, F


# ---------------------------------------------------------------- adjacency chain (DFS)
def chain_pack(container, areas, time_cap=30, verbose=True):
    """Backtracking DFS packing in which each piece touches the previous one.
    Returns the deepest chain found (list of (poly, ra))."""
    Cont = Polygon(container); best_chain = [None]; t0 = time.time()

    def candidates(prev_poly, placed_union, area):
        L = sqrt(2 * area); os = orient_variants(L)
        anchors = free_vertices(Cont if prev_poly is None else prev_poly)
        cont_p = prep(Cont)
        out = []
        for v in anchors:
            vv = np.array(v)
            for o, ra in os:
                poly = Polygon(o + vv)
                if not poly.is_valid: continue
                if not cont_p.contains(poly.buffer(-TOL)): continue
                if placed_union is not None and poly.buffer(-TOL).intersects(placed_union):
                    continue
                if prev_poly is not None and poly.distance(prev_poly) > 1e-7:
                    continue                              # must touch predecessor
                sc = (poly.boundary.intersection(
                        (prev_poly.boundary if prev_poly is not None else Cont.boundary)).length)
                out.append((sc, o + vv, ra, poly))
        out.sort(key=lambda z: -z[0])
        return out

    def dfs(i, prev_poly, placed_union, chain):
        if len(chain) > len(best_chain[0] or []):
            best_chain[0] = list(chain)
        if i >= len(areas) or (time_cap and time.time() - t0 > time_cap):
            return
        for _, pts, ra, poly in candidates(prev_poly, placed_union, areas[i]):
            u = poly if placed_union is None else placed_union.union(poly)
            chain.append((pts, ra))
            dfs(i + 1, poly, u, chain)
            chain.pop()
            if time_cap and time.time() - t0 > time_cap:
                return

    dfs(0, None, None, [])
    if verbose:
        print(f"  chain: {len(best_chain[0] or [])} pieces, {time.time()-t0:.0f}s")
    return best_chain[0] or []


# ---------------------------------------------------------------- integer dissection
def affine_unit_to(tri_pts, ra_index):
    """Affine map A(x)=M x + t sending UNIT_T (RA at vertex 0, legs to v1,v2) onto the given
    triangle (RA at ra_index). Preserves the leg order, so the child packing maps in similarly."""
    src = np.array(UNIT_T)
    order = [ra_index, (ra_index + 1) % 3, (ra_index + 2) % 3]
    dst = np.array(tri_pts)[order]
    # solve M (2x2), t (2,) from 3 correspondences
    Asys = np.zeros((6, 6)); bsys = np.zeros(6)
    for i in range(3):
        x, y = src[i]
        Asys[2*i]   = [x, y, 0, 0, 1, 0]; bsys[2*i]   = dst[i, 0]
        Asys[2*i+1] = [0, 0, x, y, 0, 1]; bsys[2*i+1] = dst[i, 1]
    sol = np.linalg.solve(Asys, bsys)
    M = np.array([[sol[0], sol[1]], [sol[2], sol[3]]]); t = np.array([sol[4], sol[5]])
    return lambda pts: (np.array(pts) @ M.T + t)

def restricted_child_areas(state_prime, P, dP, Rlt, kmax):
    """Areas of the children [state*q], q>=state_prime, RELATIVE to the parent cell.
    area[state*q]/area[state] = (1/q) prod_{state<=r<q}(1-1/r)."""
    idx0 = int(np.searchsorted(P, state_prime))
    out = []
    Rbase = Rlt[idx0]                                    # prod_{r<state}
    for j in range(idx0, min(idx0 + kmax, len(P))):
        q = int(P[j])
        rel = (Rlt[j] / Rbase) / q                       # (1/q) prod_{state<=r<q}
        out.append((q, rel))
    return out

def _tri_area(pts):
    p = np.asarray(pts); return 0.5 * abs((p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]))

def integer_dissection(P, dP, Rlt, level1=42, max_level=3, min_rel=0.028, min_area=2e-3, kids=48):
    """Recursive least-prime-factor dissection with EXACT Markov self-similarity (Prop. 'markov'):
    the sub-dissection of a cell depends only on its state s = P(cell). We pack the children of
    state s once, in the unit triangle (a fixed relative pattern, children with relative area
    >= min_rel), cache it, and map that same pattern into every cell of state s by the similarity
    carrying the unit triangle onto the cell. Thus [5], [15]=3*5, [45]=9*5, [75]=15*5, ... are all
    the identical pattern up to the scaling/rotation of their containing cells (e.g. [15] is a
    1/sqrt(3) copy of [5]). A cell is drawn subdivided (a 'shell', bracket-labelled) when it is
    large enough and below max_level; otherwise it is a leaf carrying its number."""
    cells = []
    state_cache = {}

    def canonical_children(state):
        """Fixed relative packing of state-s children in the unit triangle (cached)."""
        if state not in state_cache:
            ch = [(q, rel) for (q, rel) in restricted_child_areas(state, P, dP, Rlt, kids)
                  if rel >= min_rel]
            placed, _ = greedy_pack(UNIT_T, [rel for _, rel in ch], verbose=False)
            state_cache[state] = [(ch[i][0], placed[i][0], placed[i][1]) for i in range(len(placed))]
        return state_cache[state]

    def recurse(tri, ra, state, label, level):
        area = _tri_area(tri)
        kids_c = canonical_children(state)
        subdivide = (level < max_level) and (area >= 4 * min_area) and len(kids_c) >= 2
        poly = [[float(x), float(y)] for x, y in tri]
        if not subdivide:
            cells.append({"poly": poly, "label": int(label), "level": level, "kind": "leaf"})
            return
        cells.append({"poly": poly, "label": int(label), "level": level, "kind": "shell",
                      "state": int(state)})
        A = affine_unit_to(tri, ra)
        for q, cpts, cra in kids_c:                 # SAME pattern for every cell of this state
            recurse(A(cpts), cra, q, label * q, level + 1)

    L1, _ = greedy_pack(UNIT_T, dP[:level1], verbose=False)
    for j, (pts, ra) in enumerate(L1):
        recurse(pts, ra, int(P[j]), int(P[j]), 1)
    return cells


# ---------------------------------------------------------------- head certificate
def head_certificate(P, dP, nmax):
    """Greedy to p_last, then decompose the free region into triangles with inradius 'widths'."""
    placed, F = greedy_pack(UNIT_T, dP[:nmax], verbose=False)
    p_last = int(P[nmax - 1])
    geoms = F.geoms if isinstance(F, MultiPolygon) else [F]
    regions, widths, areas = [], [], []
    for g in geoms:
        for tri in _triangulate(g):
            poly = Polygon(tri)
            if poly.area < 1e-12: continue
            regions.append([list(map(float, pt)) for pt in tri])
            widths.append(float(2 * poly.area / poly.length))   # ~inradius proxy
            areas.append(float(poly.area))
    return {"regions": regions, "widths": widths, "region_areas": areas, "p_last": p_last}

def _triangulate(poly):
    try:
        from shapely.ops import triangulate
        return [list(t.exterior.coords)[:3] for t in triangulate(poly) if t.within(poly.buffer(1e-9))]
    except Exception:
        c = np.array(poly.exterior.coords)[:-1]
        return [[c[0], c[i], c[i+1]] for i in range(1, len(c) - 1)]


# ---------------------------------------------------------------- writers / CLI
def _pieces_json(placed, P):
    return [{"poly": pts.tolist(), "p": int(P[i]), "d": float(0.5 * _leg2(pts))}
            for i, (pts, ra) in enumerate(placed)]

def _leg2(pts):
    e = [np.linalg.norm(pts[(i+1) % 3] - pts[i]) for i in range(3)]
    L = min(e)                       # a leg (hypotenuse is the long side)
    return L * L                     # area = 1/2 leg^2

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default="all",
                    choices=["all", "greedy", "chain", "integers", "certificate"])
    ap.add_argument("--nmax", type=int, default=355)
    ap.add_argument("--chain-nmax", type=int, default=60)
    ap.add_argument("--time-cap", type=float, default=0.0)
    a = ap.parse_args()
    P = sieve(200000); dP, Rlt = lpf_densities(P)

    if a.target in ("all", "greedy"):
        placed, F = greedy_pack(UNIT_T, dP[:a.nmax], time_cap=a.time_cap)
        json.dump(_pieces_json(placed, P), open("pieces_greedy.json", "w"))
        print(f"pieces_greedy.json: {len(placed)} pieces, primes 2..{int(P[len(placed)-1])}, "
              f"covered {1-F.area:.4f}")
    if a.target in ("all", "chain"):
        chain = chain_pack(UNIT_T, dP[:a.chain_nmax], time_cap=max(a.time_cap, 30))
        json.dump(_pieces_json(chain, P), open("pieces_chain_dfs.json", "w"))
        print(f"pieces_chain_dfs.json: {len(chain)} pieces, primes 2..{int(P[len(chain)-1]) if chain else '-'}")
    if a.target in ("all", "integers"):
        cells = integer_dissection(P, dP, Rlt)
        json.dump(cells, open("integers.json", "w"))
        print(f"integers.json: {len(cells)} cells across 3 levels")
    if a.target in ("all", "certificate"):
        H = head_certificate(P, dP, a.nmax)
        json.dump(H, open("head_certificate.json", "w"))
        print(f"head_certificate.json: {len(H['regions'])} cells at p_last={H['p_last']}")

if __name__ == "__main__":
    main()
