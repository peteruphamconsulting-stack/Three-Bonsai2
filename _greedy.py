"""Greedy triangle packer — the placement core extracted verbatim from make_data.py
(P. Upham, companion code), used by triangle_base_exact.py to reproduce the base
packing whose exact validity is certified. Float geometry decides the combinatorial
placement; the exact verifier realizes the same choices in Q(sqrt2,...,sqrt47)."""
import time
from math import sqrt, pi
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely.prepared import prep
from shapely.ops import unary_union

SQ2 = sqrt(2.0)
TOL = 1e-9
UNIT_T = [(0.0, 0.0), (SQ2, 0.0), (0.0, SQ2)]          # right angle at vertex 0


def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)


def lpf_densities(P):
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0 / P.astype(float)))[:-1]])
    Rlt = np.exp(logfac)
    return Rlt / P, Rlt


def orient_variants(L):
    base = np.array([(0.0, 0.0), (L, 0.0), (0.0, L)])   # right angle at index 0
    out = []
    for k in range(8):
        th = k * pi / 4; c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s], [s, c]])
        v = base @ R.T
        for a in range(3):
            shifted = v - v[a]
            out.append((shifted, 0, k, a))              # (verts, ra_base=0, rot k, anchor a)
    uniq = []
    for o, ra, k, a in out:
        if not any(np.allclose(o, u, atol=1e-12) for u, _, _, _ in uniq):
            uniq.append((o, ra, k, a))
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


def greedy_pack(container, areas, time_cap=0, verbose=False):
    """Corner-anchored deepest-nestle greedy. Returns list of
    (poly(3x2 array), ra_index, rot_k, anchor_a) and the final free region."""
    F = Polygon(container); placed = []; cache = {}; t0 = time.time()
    for idx, area in enumerate(areas):
        L = sqrt(2 * area); key = round(L, 12)
        os = cache.get(key)
        if os is None: os = cache[key] = orient_variants(L)
        pf = prep(F); best = None; bestsc = -1.0
        for v in free_vertices(F):
            vv = np.array(v)
            for o, ra, k, a in os:
                poly = Polygon(o + vv)
                if not poly.is_valid: continue
                if pf.contains(poly.buffer(-TOL)):
                    inter = poly.boundary.intersection(F.boundary)
                    sc = inter.length if not inter.is_empty else 0.0
                    if sc > bestsc + 1e-12:
                        bestsc = sc; best = (o + vv, ra, k, a, poly)
        if best is None:
            if verbose: print(f"  greedy jam at piece {idx}")
            break
        pts, ra, k, a, poly = best
        placed.append((pts, ra, k, a))
        F = _clean(F.difference(poly))
        if time_cap and time.time() - t0 > time_cap:
            break
    return placed, F
