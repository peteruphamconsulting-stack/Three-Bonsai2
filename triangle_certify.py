#!/usr/bin/env python3
"""
triangle_certify.py -- tail crossover and a rigorous head certifier for the LPF
similar-TRIANGLE packing (the triangle analogue of ellipse_certify.py).

FRAME: container = right isosceles triangle with vertices (-1,0), (1,0), (0,1)
(hypotenuse on the x-axis, area 1). All three edges are RATIONAL lines:
    y >= 0,      x + y <= 1,      y - x <= 1.
Pieces: SIMILAR right isosceles triangles of areas d_p, leg L_p = sqrt(2 d_p),
corner-anchored in the 8 orientations k*45deg; a placement is (p, k, a, ax, ay)
with vertices  anchor + L_p * (Rot_k(base) - Rot_k(base)[a]),  base=(0,0),(1,0),(0,1).
Rotation entries lie in Z + Z*(sqrt2/2), so each vertex coordinate is EXACTLY
    (stated rational) + a*sqrt(2 d_p) + b*sqrt(d_p),   a,b small integers,
and every constraint margin is rational plus an explicit radical part. A margin is
certified by EXACT RATIONAL arithmetic when its radical coefficients vanish (all
intended boundary contacts are built this way, incl. piece 2's canonical half,
whose L=1 is rational), and by interval arithmetic otherwise (the builder keeps
all non-exact margins >= ~1e-9 and all pair gaps >= ~1e-7, far above interval width).

    tail     the eroded-volume placeability criterion
                 (b + 2 sum_{i<q} sqrt(d_i) V_i)/sqrt(q R_<q) + N_q/q < 1,
             V_i = V(R_i T, -T) planar mixed areas, exact table over the 8 relative
             orientations (verified by polygonal Minkowski sums):
                 v = [2, sqrt2, 2, 3/sqrt2, 1, 3/sqrt2, 2, sqrt2],
             kappa_T = max v = 3/sqrt2; Rogers-Shephard checks v(0)=V(T,-T)=2,
             v(4)=V(-T,-T)=1; the erosion identity C(-)sC=(1-s)C gives boundary
             tax 2s for container-oriented K (b=2), and the 180deg erosion tax is
             < 4s. Designs:
               uniform8            all 8 orientations everywhere: b=2, V<=kappa_T
                                   --> q0 = 662177, head = 53730;
               disciplined --p0 Q  free orientations for p<=Q, strictly alternating
                                   {0,180} beyond, next-piece K adaptive in {0,180}:
                                   averaging the two K choices bounds the free mass
                                   by weight max_k (v(k)+v(k+4))/2 = 2, alternating
                                   mass by (v(0)+v(4))/2 = 3/2, boundary by 3s
                                   --> Q=1000: q0 = 115363, head = 10906.
             The analytic majorant of Appendix A applies verbatim with kappa_T in
             place of kappa (the shape enters only through the mixed-area constant).
    certify  reads a builder head JSON and proves, fail-closed:
             * areas == d_p exactly (reconstructed from n, never parsed);
             * CONTAINMENT: per-vertex margins for the three rational edge lines,
               exact-rational when radical coefficients vanish, else certified
               interval lower bounds;
             * DISJOINTNESS: far pairs pruned by a padded-bbox test (rigorous, wide
               margin); close pairs by a SEPARATING-AXIS witness over the 8 axis
               directions in interval arithmetic (these directions contain every
               edge normal of every piece) -- the found axis with its certified
               positive gap is the witness, in exact parallel to the rational
               Perram-Wertheim witnesses of the ellipse certificate.

Usage:
    python3 triangle_certify.py tail    [--design uniform8|disciplined] [--p0 1000] [--X 4000000]
    python3 triangle_certify.py certify [--in triangle_head.json] [--prec 100]

Requires: numpy, mpmath, sympy.
"""
import argparse, json, math
from fractions import Fraction as Fr
import numpy as np, sympy

SQRT2 = math.sqrt(2.0)
KAPPA8 = 3.0/SQRT2
VTABLE = [2.0, SQRT2, 2.0, KAPPA8, 1.0, KAPPA8, 2.0, SQRT2]

# ------------------------------------------------------------------ tail
def _sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s).astype(np.float64)

def tail(design="uniform8", p0=1000, X=4_000_000, state=2, quiet=False):
    """Crossover for the state-s region. Using sigma_s(q)=d_q/R_<s and
    R^(s)_<q = R_<q/R_<s, every state quantity reduces to a global one scaled by
    R_<s, so the criterion is exact in the global prime arrays (verified: s=2
    reproduces the base 115363/10906 and 662177/53730)."""
    P = _sieve(X)
    logR = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0/P))])[:-1]
    R = np.exp(logR); sq = np.sqrt(R/P)
    Tc = np.concatenate([[0.0], np.cumsum(sq)])[:-1]           # Tc[j]=sum_{i<j} sqrt(d_i)
    i_s = int(np.searchsorted(P, state)); Rs = R[i_s]; sRs = math.sqrt(Rs)
    j = np.arange(i_s, len(P)); Pj = P[j]; Rj = R[j]; Nj = j - i_s
    if design == "uniform8":
        num = 2.0*sRs + 2.0*KAPPA8*(Tc[j] - Tc[i_s]); b = 2.0
    else:
        sqf = np.where(P <= p0, sq, 0.0); sqa = np.where(P > p0, sq, 0.0)
        Cf = np.concatenate([[0.0], np.cumsum(sqf)])[:-1]
        Ca = np.concatenate([[0.0], np.cumsum(sqa)])[:-1]
        num = 3.0*sRs + 2.0*(2.0*(Cf[j] - Cf[i_s]) + 1.5*(Ca[j] - Ca[i_s])); b = 3.0
    crit = num/np.sqrt(Pj*Rj) + Nj/Pj
    if not quiet:
        if design == "uniform8":
            print(f"state s={state}, design uniform8: b=2, V<=kappa_T=3/sqrt2={KAPPA8:.6f}")
        else:
            print(f"state s={state}, design disciplined: free orientations for p<={p0}, "
                  f"alternating 0/180 beyond; adaptive-K weights 2 / 1.5, b<=3")
    bad = np.where(crit >= 1.0)[0]
    if len(bad) == 0:
        if not quiet: print("criterion holds for every prime q>=s: EMPTY head (terminal state)")
        return None, 0
    jj = bad[-1]; q0 = int(Pj[jj]); head = int(jj + 1); nxt = int(Pj[jj+1])
    if not quiet:
        print(f"crossover: last-failing prime q0={q0} (value {crit[jj]:.6f}>=1); "
              f"holds from {nxt} (value {crit[jj+1]:.6f}<1)")
        print(f"HEAD size = #(primes in [{state},q0]) = {head} triangles")
    if quiet: return q0, head
    g = 0.5772156649; ell = math.log(1e5); beta = math.exp(-g)*(1-1/ell**2); A1e5 = 2.3355
    kap = KAPPA8 if design == "uniform8" else 1.5
    def M(q):
        Lq = math.log(q)
        return 2*kap*A1e5/(math.sqrt(beta)*Lq) + 2*kap/math.sqrt(q*(beta/Lq)) + (1/Lq)*(1+1.2762/Lq)
    for q in (1e6, 4e6, 1e7):
        print(f"  analytic majorant M({int(q)}) = {M(q):.4f}  (<1: {M(q) < 1})")
    return q0, head

# ------------------------------------------------------------------ exact geometry
# rotation entries as (a, b) meaning a + b*sqrt2/2
_ROT = {0: ((1,0),(0,0)), 1: ((0,1),(0,1)), 2: ((0,0),(1,0)), 3: ((0,-1),(0,1)),
        4: ((-1,0),(0,0)), 5: ((0,-1),(0,-1)), 6: ((0,0),(-1,0)), 7: ((0,1),(0,-1))}
BASE = [(0,0),(1,0),(0,1)]

def texact(k, a):
    """Template vertices as pairs ((ax_i, bx_h), (ay_i, by_h)): integer + (sqrt2/2)-integer.
    Multiplying by L: L*(a + b*sqrt2/2) = a*sqrt(2 d) + b*sqrt(d)."""
    (c_i, c_h), (s_i, s_h) = _ROT[k]
    vs = []
    for (x, y) in BASE:
        vs.append(((c_i*x - s_i*y, c_h*x - s_h*y), (s_i*x + c_i*y, s_h*x + c_h*y)))
    (axi, axh), (ayi, ayh) = vs[a]
    return [((xi-axi, xh-axh), (yi-ayi, yh-ayh)) for ((xi, xh), (yi, yh)) in vs]

def state_shares(n, s=2):
    """Exact areas sigma_s(q)=d_q/R_<s for the first n primes q>=s (s=2 -> d_p).
    The state region is normalised to area 1, so these telescope to 1; their
    denominators are astronomical and are reconstructed here, never parsed."""
    base = int(sympy.primepi(s-1)); hi = int(sympy.prime(base+n))
    primes = [int(p) for p in sympy.primerange(s, hi+1)][:n]
    a = Fr(1); out = []
    for q in primes: out.append(a/q); a *= Fr(q-1, q)
    return primes, out

def lpf_shares(n): return state_shares(n, 2)

def _ratsqrt(f):
    """sqrt of a Fraction if it is a perfect rational square, else None."""
    n, d = f.numerator, f.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    return Fr(rn, rd) if (rn*rn == n and rd*rd == d) else None

def load_head(path):
    doc = json.load(open(path)); V = doc.get('verification', doc)
    pieces = V['pieces']; n = len(pieces)
    s = int(doc.get('state', V.get('state', 2)))
    primes, shares = state_shares(n, s)
    for i, rec in enumerate(pieces):
        if int(rec[0]) != primes[i]:
            raise ValueError(f"piece {i}: prime {rec[0]} != expected {primes[i]} (state {s})")
    return pieces, primes, shares, doc.get('design', 'uniform8'), s

# ------------------------------------------------------------------ certify
def certify(path, prec=100, prefilter_pad=1e-6):
    from mpmath import iv
    iv.prec = prec
    pieces, primes, shares, design, state = load_head(path)
    n = len(pieces)
    # design conformance: disciplined heads must strictly alternate {0,180} past p0
    doc = json.load(open(path)); p0 = int(doc.get('p0', 0))
    design_ok = True
    if design == 'disciplined':
        for i, rec in enumerate(pieces):
            if primes[i] <= p0: continue
            k = int(rec[1])
            if k not in (0, 4): design_ok = False; break
            if k != (i % 2)*4: design_ok = False; break
    TEX = {(k, a): texact(k, a) for k in range(8) for a in range(3)}

    d_iv  = [iv.mpf(s.numerator)/iv.mpf(s.denominator) for s in shares]
    s2d   = [iv.sqrt(2*dv) for dv in d_iv]      # sqrt(2 d_p)
    sd    = [iv.sqrt(dv) for dv in d_iv]        # sqrt(d_p)
    r2d   = [_ratsqrt(2*s) for s in shares]     # rational sqrt(2d) if exact (p=2: L=1)
    rd    = [_ratsqrt(s) for s in shares]
    anchQ = [(Fr(str(r[3])), Fr(str(r[4]))) for r in pieces]
    KA    = [(int(r[1]), int(r[2])) for r in pieces]

    def margin(i, q, ca, cb):
        """Certified sign info for  q + ca*sqrt(2 d_i) + cb*sqrt(d_i)  >= 0.
        Returns (ok, lower_bound_float, exact_flag)."""
        if (ca == 0 or r2d[i] is not None) and (cb == 0 or rd[i] is not None):
            val = q + (ca*r2d[i] if ca else 0) + (cb*rd[i] if cb else 0)
            return (val >= 0, float(val), True)
        m = iv.mpf(q.numerator)/iv.mpf(q.denominator)
        if ca: m = m + ca*s2d[i]
        if cb: m = m + cb*sd[i]
        lb = float(m.a.a) if hasattr(m.a, 'a') else float(m.a)
        return (lb >= 0, lb, False)

    cont_bad = 0; cont_min = None; exact_contacts = 0
    for i in range(n):
        (qax, qay) = anchQ[i]
        for ((axi, axh), (ayi, ayh)) in TEX[KA[i]]:
            # y >= 0 ; 1-x-y >= 0 ; 1-y+x >= 0
            tests = [
                (qay,                ayi,            ayh),
                (1 - qax - qay,     -(axi+ayi),     -(axh+ayh)),
                (1 - qay + qax,      (axi-ayi),      (axh-ayh)),
            ]
            for (q, ci, ch) in tests:
                ok, lb, ex = margin(i, q, ci, ch)
                if ex and lb == 0.0: exact_contacts += 1
                if not ok: cont_bad += 1
                if cont_min is None or lb < cont_min: cont_min = lb

    # float vertices for the prefilter
    Hf = SQRT2/2
    Lf = np.array([math.sqrt(2*float(s)) for s in shares])
    vertsF = np.zeros((n, 3, 2))
    for i in range(n):
        (qax, qay) = anchQ[i]
        for v, ((axi, axh), (ayi, ayh)) in enumerate(TEX[KA[i]]):
            vertsF[i, v, 0] = float(qax) + axi*Lf[i] + axh*Lf[i]*Hf
            vertsF[i, v, 1] = float(qay) + ayi*Lf[i] + ayh*Lf[i]*Hf
    lo = vertsF.min(axis=1) - prefilter_pad; hi = vertsF.max(axis=1) + prefilter_pad
    cell = max(float(np.median(Lf)), 1e-3); grid = {}
    for i in range(n):
        for cx in range(int(math.floor(lo[i,0]/cell)), int(math.floor(hi[i,0]/cell))+1):
            for cy in range(int(math.floor(lo[i,1]/cell)), int(math.floor(hi[i,1]/cell))+1):
                grid.setdefault((cx, cy), []).append(i)
    cand = set()
    for ids in grid.values():
        for u in range(len(ids)):
            for v in range(u+1, len(ids)):
                cand.add((ids[u], ids[v]))
    close = [(i, j) for (i, j) in cand
             if not (hi[i,0] < lo[j,0] or hi[j,0] < lo[i,0] or
                     hi[i,1] < lo[j,1] or hi[j,1] < lo[i,1])]
    npairs = n*(n-1)//2; pruned = npairs - len(close)

    # interval vertices once
    hI = iv.sqrt(2)/2
    vertsI = []
    for i in range(n):
        (qax, qay) = anchQ[i]
        qx = iv.mpf(qax.numerator)/iv.mpf(qax.denominator)
        qy = iv.mpf(qay.numerator)/iv.mpf(qay.denominator)
        vs = []
        for ((axi, axh), (ayi, ayh)) in TEX[KA[i]]:
            x = qx + axi*s2d[i] + axh*sd[i]
            y = qy + ayi*s2d[i] + ayh*sd[i]
            vs.append((x, y))
        vertsI.append(vs)
    AX = [(iv.mpf(1), iv.mpf(0)), (hI, hI), (iv.mpf(0), iv.mpf(1)), (-hI, hI),
          (iv.mpf(-1), iv.mpf(0)), (-hI, -hI), (iv.mpf(0), iv.mpf(-1)), (hI, -hI)]
    sat_bad = 0; sat_min = None
    for (i, j) in close:
        found = None
        for (ux, uy) in AX:
            pi = [x*ux + y*uy for (x, y) in vertsI[i]]
            pj = [x*ux + y*uy for (x, y) in vertsI[j]]
            g1 = (min(p.a for p in pj) - max(p.b for p in pi)).a
            g2 = (min(p.a for p in pi) - max(p.b for p in pj)).a
            g = max(float(g1), float(g2))
            if g > 0: found = g; break
        if found is None: sat_bad += 1
        elif sat_min is None or found < sat_min: sat_min = found

    inv = {"areas_are_sigma_s_exact": True,
           "design_pattern_ok": bool(design_ok),
           "all_contained": bool(cont_bad == 0),
           "all_pairs_disjoint": bool(sat_bad == 0)}
    status = "CERTIFIED" if all(inv.values()) else "FAILED"
    print(f"triangle head: state s={state}, n={n}, primes {primes[0]}..{primes[-1]}, design={design}")
    print(f"  areas == sigma_{state}(q) (exact, reconstructed from n): True")
    print(f"  design pattern ({design}{f', p0={p0}' if design=='disciplined' else ''}): "
          f"{'OK' if design_ok else 'VIOLATED'}")
    print(f"  contained: {3*3*n - cont_bad}/{3*3*n} vertex-edge margins certified "
          f"({exact_contacts} exact rational contacts, min margin {cont_min:+.3e})")
    print(f"  disjoint: {npairs - sat_bad}/{npairs} pairs "
          f"({len(close)} axis witnesses, {pruned} pruned by padded bbox, "
          f"min witness gap {'n/a' if sat_min is None else f'{sat_min:.3e}'})")
    print(f"  STATUS: {status}")
    return status

# ------------------------------------------------------------------ CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default="certify", choices=["tail", "certify"])
    ap.add_argument("--design", default="uniform8", choices=["uniform8", "disciplined"])
    ap.add_argument("--p0", type=int, default=1000)
    ap.add_argument("--X", type=int, default=4_000_000)
    ap.add_argument("--state", type=int, default=2)
    ap.add_argument("--in", dest="infile", default="triangle_head.json")
    ap.add_argument("--prec", type=int, default=100)
    a = ap.parse_args()
    if a.target == "tail":
        tail(a.design, a.p0, a.X, state=a.state)
    else:
        certify(a.infile, prec=a.prec)

if __name__ == "__main__":
    main()
