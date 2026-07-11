#!/usr/bin/env python3
"""
triangle_base_exact.py  --  EXACT verifier for the base packing of Proposition 6.3
("A Geometric Sieve of Eratosthenes", P. Upham).

Certifies, for the greedy base packing of the primes 2,3,...,47 (15 similar right
isosceles triangles of areas d_2,...,d_47 inside the container T = (0,0),(sqrt2,0),
(0,sqrt2)):

  * every coordinate lies in Q(sqrt2, sqrt3, ..., sqrt47), reconstructed exactly from
    the greedy's combinatorial choices (anchor vertex + 45-deg rotation + edge length);
  * each piece has area EXACTLY d_p                          (15 exact area identities);
  * every piece lies inside T                                (135 = 15*9 exact sign-tests);
  * the 15 pieces have pairwise-disjoint interiors           (105 = C(15,2) exact
    separating-edge certificates);
  * the free region has area EXACTLY R_{<50} = prod_{p<50}(1-1/p), and it CONTAINS a disk
    of certified radius >= 0.106 -- a rigorous lower bound on its inradius. The disk centre
    is proposed by GEOS in floating point, then certified to lie in the exact free region,
    and its radius is a floored mpmath-interval lower bound on the centre's exact distance
    to the boundary (we do not claim GEOS found the true maximiser).

Float geometry (shapely) is used only to fix the combinatorics (which anchor, which
orientation, which triangulation); every certified inequality is then decided exactly
in Q(sqrt2,...,sqrt47) (sympy) or by directed-rounding interval arithmetic (mpmath).

Requires: numpy, sympy, mpmath, shapely.  Deterministic; ~a few seconds.
Usage:  python3 triangle_base_exact.py
"""
import numpy as np, sympy as sp
from mpmath import iv
import _greedy as G

S2 = sp.sqrt(2)
DIRS = [(sp.Integer(1), sp.Integer(0)), (sp.Integer(-1), sp.Integer(0)),
        (sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(-1)),
        ( S2/2,  S2/2), ( S2/2, -S2/2), (-S2/2,  S2/2), (-S2/2, -S2/2)]
DIRS_F = [(float(x), float(y)) for x, y in DIRS]


def esign(e):
    """Exact sign (-1/0/+1) of a Q(sqrt2,...,sqrt47) element. Fails CLOSED: it never
    substitutes floating-point evidence for a sign. For a nonzero element a directed-
    rounding mpmath interval at rising precision resolves the sign rigorously; a genuine
    zero is caught symbolically. If neither resolves it, we raise rather than guess."""
    e = sp.simplify(e)
    if e.is_zero:      return 0
    if e.is_positive:  return 1
    if e.is_negative:  return -1
    if sp.simplify(sp.radsimp(e)) == 0:                   # second symbolic zero test
        return 0
    for pr in (200, 400, 800, 1600, 3200):               # rigorous interval enclosure
        try:
            x = to_iv(e, pr)
        except Exception:
            break
        if x.a > 0: return 1
        if x.b < 0: return -1
    raise ArithmeticError(f"esign: could not certify the sign of {e!r}")


def _key(pt, nd=8):
    return (round(float(pt[0]), nd), round(float(pt[1]), nd))


def reconstruct():
    """Return P, exact triangles [(p, [(x,y)*3], d_p)], the greedy placements, free
    region F (shapely), the exact vertex registry, and a frozen combinatorial record
    (per piece: prime, anchor index, rotation, per-edge direction+magnitude, and the
    resulting exact vertices as sympy srepr) that reconstruct_from_record() can re-verify
    with no float greedy / GEOS dependence."""
    P = G.sieve(200); dP, _ = G.lpf_densities(P)
    placed, F = G.greedy_pack(G.UNIT_T, dP[:15])
    reg = {_key((0, 0)): (sp.Integer(0), sp.Integer(0)),
           _key((G.SQ2, 0)): (S2, sp.Integer(0)),
           _key((0, G.SQ2)): (sp.Integer(0), S2)}
    R = sp.Integer(1); tris = []; rec = []
    for i, (pts, ra, k, a) in enumerate(placed):
        p = int(P[i])
        d = R / p                                     # d_p = R_{<p}/p  (R already prod_{q<p})
        R = R * sp.Rational(p - 1, p)                 # advance to R_{<next}
        Lp = sp.sqrt(2 * d); Lpf = float(Lp)
        anchor = reg[_key(pts[a])]
        verts = [None, None, None]; verts[a] = anchor
        edges = []
        for j in range(3):
            if j == a: continue
            e = np.array(pts[j]) - np.array(pts[a]); elen = float(np.hypot(*e))
            is_hyp = abs(elen - Lpf) >= 1e-6
            mag = Lp * S2 if is_hyp else Lp
            u = e / elen
            di = min(range(8), key=lambda t: (u[0]-DIRS_F[t][0])**2 + (u[1]-DIRS_F[t][1])**2)
            dx, dy = DIRS[di]
            verts[j] = (sp.simplify(anchor[0] + mag*dx), sp.simplify(anchor[1] + mag*dy))
            reg[_key(pts[j])] = verts[j]
            edges.append({'vertex': j, 'dir_index': int(di), 'hypotenuse': bool(is_hyp)})
        tris.append((p, verts, d))
        rec.append({'p': p, 'anchor': int(a), 'rotation': int(k), 'edges': edges,
                    'verts': [[sp.srepr(vx), sp.srepr(vy)] for (vx, vy) in verts]})
    return P, tris, placed, F, reg, rec


def reconstruct_from_record(rec):
    """Rebuild the exact base packing from a FROZEN record (no float greedy, no GEOS for the
    combinatorics). The record's exact vertices (sympy srepr) are the authoritative data; all
    exact/interval certificates are then re-run on them. GEOS is used only to build the free
    region F for the disk-centre proposal, which is certified independently afterward."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    P = [e['p'] for e in rec]
    reg = {_key((0, 0)): (sp.Integer(0), sp.Integer(0)),
           _key((G.SQ2, 0)): (S2, sp.Integer(0)),
           _key((0, G.SQ2)): (sp.Integer(0), S2)}
    R = sp.Integer(1); tris = []; polys = []
    for e in rec:
        p = int(e['p']); d = R / p; R = R * sp.Rational(p - 1, p)
        verts = [(sp.sympify(vx), sp.sympify(vy)) for (vx, vy) in e['verts']]
        for v in verts:
            reg[_key((float(v[0]), float(v[1])))] = v
        tris.append((p, verts, d))
        polys.append(Polygon([(float(x), float(y)) for (x, y) in verts]))
    cont = Polygon([(0.0, 0.0), (G.SQ2, 0.0), (0.0, G.SQ2)])
    F = cont.difference(unary_union(polys))
    return P, tris, None, F, reg


# ---------------------------------------------------------------- exact certificates
def check_areas(tris):
    fails = 0
    for p, v, d in tris:
        (x0, y0), (x1, y1), (x2, y2) = v
        A = sp.Rational(1, 2) * ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0))
        if sp.simplify(sp.Abs(A) - d) != 0:
            fails += 1; print(f"  area mismatch p={p}")
    return fails == 0


def check_containment(tris):
    """135 sign-tests: each of 15 pieces x 3 vertices x 3 container half-planes
    (x>=0, y>=0, sqrt2 - x - y >= 0)."""
    n = 0; passed = 0
    for p, v, d in tris:
        for (x, y) in v:
            for expr in (x, y, S2 - x - y):
                n += 1
                if esign(expr) >= 0: passed += 1
                else: print(f"  containment FAIL p={p} at ({x},{y})")
    return n, passed


def check_disjoint(tris):
    """C(15,2)=105 separating-edge certificates: for each pair, exhibit an edge of one
    triangle whose supporting line has that triangle's own interior strictly on one side
    (fixed by its third vertex) and the *entire other triangle* on the opposite closed
    side. That is the actual separating-axis condition for disjoint interiors."""
    def cross(Va, Vb, px, py):
        (ax, ay), (bx, by) = Va, Vb
        return (bx-ax)*(py-ay) - (by-ay)*(px-ax)
    def separates(Va, Vb, Vc, other):
        # edge Va->Vb of a triangle whose third (interior-side) vertex is Vc.
        s_own = esign(cross(Va, Vb, Vc[0], Vc[1]))
        if s_own == 0:
            return False                                   # degenerate: line passes through Vc
        for (px, py) in other:                             # every other vertex on the far closed side
            if s_own * esign(cross(Va, Vb, px, py)) > 0:   # strictly on the own-interior side -> not separating
                return False
        return True
    pairs = 0; certified = 0
    for i in range(len(tris)):
        for j in range(i+1, len(tris)):
            pairs += 1
            Vi = tris[i][1]; Vj = tris[j][1]
            ok = False
            # candidate separators: the 3 edges of each triangle, each tagged with its opposite vertex
            cands = ([(Vi[a], Vi[(a+1) % 3], Vi[(a+2) % 3], Vj) for a in range(3)] +
                     [(Vj[a], Vj[(a+1) % 3], Vj[(a+2) % 3], Vi) for a in range(3)])
            for (Va, Vb, Vc, other) in cands:
                if separates(Va, Vb, Vc, other):
                    ok = True; break
            if ok: certified += 1
            else: print(f"  no separating edge for pieces {tris[i][0]},{tris[j][0]}")
    return pairs, certified


# ---------------------------------------------------------------- free region widths
def to_iv(expr, prec=200):
    """Rigorous mpmath-interval enclosure of a Q(sqrt2,...) element (recursive)."""
    iv.prec = prec
    expr = sp.sympify(expr)
    if expr.is_Integer:  return iv.mpf(int(expr))
    if expr.is_Rational: return iv.mpf(int(expr.p)) / iv.mpf(int(expr.q))
    if isinstance(expr, sp.Pow):
        if expr.exp == sp.Rational(1, 2): return iv.sqrt(to_iv(expr.base, prec))
        if expr.exp.is_Integer:           return to_iv(expr.base, prec) ** int(expr.exp)
    if expr.is_Add:
        r = iv.mpf(0)
        for a in expr.args: r = r + to_iv(a, prec)
        return r
    if expr.is_Mul:
        r = iv.mpf(1)
        for a in expr.args: r = r * to_iv(a, prec)
        return r
    raise ValueError(f"cannot enclose {expr!r}")


def _exact_vertex(reg, ptf):
    return reg[_key(ptf)]


def _tri_inradius_iv(v, prec=200):
    """Rigorous mpmath-interval inradius (= area/semiperimeter) of a triangle with
    exact-sympy vertices v."""
    pts = [(to_iv(x, prec), to_iv(y, prec)) for (x, y) in v]
    def dist(p, q): return iv.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)
    A = abs((pts[1][0]-pts[0][0])*(pts[2][1]-pts[0][1]) -
            (pts[2][0]-pts[0][0])*(pts[1][1]-pts[0][1])) / 2
    s = (dist(pts[0], pts[1]) + dist(pts[1], pts[2]) + dist(pts[2], pts[0])) / 2
    return A / s


def free_region(P, tris, F, reg, prec=200):
    """Exact free area (= R_{<50}) and a rigorously certified inscribed-disk radius."""
    # exact identity: free area = 1 - sum_{p<=47} d_p = R_{<53} = R_{<50}
    dsum = sp.Integer(0)
    for _, _, d in tris: dsum += d
    free_exact = sp.nsimplify(1 - dsum)
    R50 = sp.Integer(1)
    for q in [int(p) for p in P if p < 50]: R50 *= sp.Rational(q - 1, q)
    area_ok = sp.simplify(free_exact - R50) == 0
    # largest inscribed disk: GEOS proposes the centre in floating point; we then (i) certify
    # rigorously that this centre lies in the exact free region, and (ii) return a floored
    # interval lower bound on its distance to the exact boundary. This certifies that F
    # CONTAINS a disk of that radius -- i.e. a lower bound on inradius(F) -- not that GEOS
    # found the true maximiser.
    import shapely
    try:
        mic = shapely.maximum_inscribed_circle(F)
        cx, cy = mic.coords[0]
    except Exception:
        from shapely import ops as shops
        c = shops.polylabel(F, tolerance=1e-7); cx, cy = c.x, c.y
    center_ok = _point_in_free(cx, cy, tris, prec)
    Rdisk = _certify_disk(cx, cy, reg, F, prec) if center_ok else 0.0
    return area_ok, float(free_exact), Rdisk, center_ok


def _seg_dist_iv(C, A, B):
    """Rigorous mpmath-interval LOWER bound on dist(C, segment AB) (A,B,C are iv pairs).
    Key fact: dist(C, seg) >= dist(C, line AB) always, since seg subset line. So the line
    distance is a valid lower bound whenever the foot might lie in [0,1]; it is sharpened
    to the exact endpoint distance only when the projection is CERTIFIED beyond an endpoint.
    Never selects a candidate by comparing upper endpoints."""
    ABx, ABy = B[0]-A[0], B[1]-A[1]
    ACx, ACy = C[0]-A[0], C[1]-A[1]
    L2 = ABx**2 + ABy**2
    t = (ACx*ABx + ACy*ABy) / L2
    dline = abs(ABx*ACy - ABy*ACx) / iv.sqrt(L2)          # distance to the infinite line
    if t.b <= 0:                                          # foot certified at/before A
        return iv.sqrt(ACx**2 + ACy**2)                  #   -> exact distance to A
    if t.a >= 1:                                          # foot certified at/after B
        return iv.sqrt((C[0]-B[0])**2 + (C[1]-B[1])**2)  #   -> exact distance to B
    return dline                                          # foot possibly interior -> line dist (rigorous lower bound)


def _certify_disk(cx, cy, reg, F, prec=200):
    """Rigorous LOWER bound on dist((cx,cy), boundary(F)): the min over boundary SEGMENTS
    of the certified lower bound _seg_dist_iv(...).a. Aggregating by the segments' LOWER
    endpoints is what makes the return value a rigorous lower bound; selecting the interval
    of smallest UPPER endpoint (as before) does not. The caller must separately certify
    (via _point_in_free) that the centre lies in F -- a point outside F can also have
    positive boundary distance, so distance alone does not certify an inscribed disk."""
    iv.prec = prec
    C = (iv.mpf(float(cx)), iv.mpf(float(cy)))
    def edges(ring):
        cc = list(ring.coords)
        return [(_exact_vertex(reg, cc[i]), _exact_vertex(reg, cc[i+1])) for i in range(len(cc)-1)]
    E = edges(F.exterior)
    for r in F.interiors: E += edges(r)
    best_lo = None
    for (Va, Vb) in E:
        A = (to_iv(Va[0], prec), to_iv(Va[1], prec))
        B = (to_iv(Vb[0], prec), to_iv(Vb[1], prec))
        lo = _seg_dist_iv(C, A, B).a                       # rigorous lower bound for this edge
        best_lo = lo if best_lo is None else (lo if lo < best_lo else best_lo)
    return float(best_lo) if best_lo is not None else 0.0


def _point_in_free(cx, cy, tris, prec=200):
    """Rigorously certify that the float point (cx,cy) lies in the INTERIOR of the free
    region: strictly inside the container T = (0,0),(sqrt2,0),(0,sqrt2) and strictly
    OUTSIDE every packed triangle. The container test uses a directed-rounding interval
    enclosure of the point; the outside-each-triangle test exhibits, per triangle, one
    edge whose exact interior side (fixed by the third vertex) has the enclosed point
    strictly on the far side. Returns True only when every test resolves with the correct
    strict sign; without it a positive boundary distance would not certify an inscribed disk."""
    iv.prec = prec
    Cx, Cy = iv.mpf(float(cx)), iv.mpf(float(cy))
    R2 = iv.sqrt(iv.mpf(2))
    if not (Cx.a > 0 and Cy.a > 0 and (R2 - Cx - Cy).a > 0):    # strictly inside container
        return False
    for (_p, v, _d) in tris:                                    # strictly outside each triangle
        outside = False
        for a in range(3):
            Va, Vb, Vc = v[a], v[(a+1) % 3], v[(a+2) % 3]
            s_own = esign((Vb[0]-Va[0])*(Vc[1]-Va[1]) - (Vb[1]-Va[1])*(Vc[0]-Va[0]))
            if s_own == 0:
                continue
            ax, ay = to_iv(Va[0], prec), to_iv(Va[1], prec)
            bx, by = to_iv(Vb[0], prec), to_iv(Vb[1], prec)
            crC = (bx-ax)*(Cy-ay) - (by-ay)*(Cx-ax)            # interval cross-product for C
            if (s_own > 0 and crC.b < 0) or (s_own < 0 and crC.a > 0):
                outside = True; break                          # C strictly on the far side of this edge
        if not outside:
            return False
    return True


# ---------------------------------------------------------------- driver
def main():
    import argparse, json, os
    ap = argparse.ArgumentParser(description="Exact Proposition 6.3 base-packing verifier.")
    ap.add_argument('--frozen', nargs='?', const='base_combinatorics.json', default=None,
                    help='verify from a frozen combinatorial record (no float greedy/GEOS for '
                         'the combinatorics); default file base_combinatorics.json')
    ap.add_argument('--emit-record', dest='emit', default=None,
                    help='run the float greedy once and write the frozen combinatorial record here')
    args = ap.parse_args()

    if args.emit:
        _, _, _, _, _, rec = reconstruct()
        json.dump(rec, open(args.emit, 'w'), indent=1)
        print(f"wrote frozen base combinatorics: {args.emit}  ({len(rec)} pieces)")
        return

    if args.frozen:
        rec = json.load(open(args.frozen))
        P, tris, placed, F, reg = reconstruct_from_record(rec)
        print(f"base packing: reconstructed {len(tris)} triangles from {args.frozen} "
              f"(frozen record; no greedy/GEOS for combinatorics)")
    else:
        P, tris, placed, F, reg, rec = reconstruct()
        print(f"base packing: {len(tris)} triangles for primes {[t[0] for t in tris]} "
              f"(float greedy; use --emit-record to freeze, --frozen to replay)")

    areas_ok = check_areas(tris)
    print(f"[exact] areas == d_p                : {areas_ok}")

    n_c, pass_c = check_containment(tris)
    print(f"[exact] containment sign-tests      : {pass_c}/{n_c} pass  (15 x 3 verts x 3 edges)")

    n_p, cert_p = check_disjoint(tris)
    print(f"[exact] separating-edge certificates: {cert_p}/{n_p} pairs  (C(15,2))")

    cells_ok, free_area, Rdisk, center_ok = free_region(P, tris, F, reg)
    print(f"[exact] free area == R_<50          : {cells_ok}   (= {free_area:.10f})")
    print(f"[interval] inscribed disk           : centre in free region = {center_ok}, "
          f"radius >= {Rdisk:.6f}  (rigorous lower bound on inradius)")

    disk_ok = bool(center_ok) and (Rdisk >= 0.106)
    ok = areas_ok and (pass_c == n_c) and (cert_p == n_p) and cells_ok and disk_ok
    print("\nCERTIFIED - base packing exactly valid in Q(sqrt2,...,sqrt47); free region contains"
          "\n           a disk of certified radius >= 0.106 (inradius >= 0.106)" if ok else "\nFAILED")


if __name__ == "__main__":
    main()
