#!/usr/bin/env python3
"""
triangle_build.py -- head builder for the LPF similar-triangle packing (discovery
tool; all rigor lives in triangle_certify.py).

FRAME: the container is the right isosceles triangle with vertices (-1,0), (1,0),
(0,1): hypotenuse (length 2) on the x-axis, right angle at the apex (0,1), area 1.
All three edges are RATIONAL lines --  y>=0,  x+y<=1,  y-x<=1  -- so a boundary
contact certifies exactly whenever the contacting vertex's radical coefficient
vanishes (template entry exactly 0 for y=0; entries summing exactly to 0 for the
legs), given a dyadic anchor. Anchors are snapped to the 2^-40 grid, making the
rational parts of all constraint margins exact in double precision.

Piece for the i-th prime p: similar right isosceles triangle of area d_p, leg
L=sqrt(2 d_p), corner-anchored in 8 orientations k*45deg x 3 anchor corners
(24 templates). Piece 2 (L=1 exactly) is seeded as the canonical half
(k=0, a=0, anchor (0,0)): the self-similar split, exact by rational arithmetic.

Method: anchors = container corners + all placed vertices, scanned freshest-first
then by priority; candidates are scored by CONTACT LENGTH (shared boundary with the
container and with neighbours -- the paper's deepest-nestle rule); overlap is a
vectorized separating-axis test over the 8 axis directions with strict gap eps;
uniform spatial grid; escalating back-off jam repair; checkpoint/resume.

Designs (must match triangle_certify.py tail):
    uniform8              all 24 templates everywhere        (head 53730)
    disciplined --p0 Q    free for p<=Q, then alternating strictly between
                          orientations 0 and 180 (head 10906 at Q=1000)

Usage:
    python3 triangle_build.py --target-n 500 --out triangle_head.json
    python3 triangle_build.py --target-n 53730 --out triangle_head_n53730.json \
        --checkpoint tri_ckpt.json --ckpt-every 500 [--resume tri_ckpt.json]

Requires: numpy.
"""
import argparse, json, math, time
import numpy as np

H = math.sqrt(2.0)/2.0
SNAP = 2.0**-40
ETA = 1e-7
COS = [1, H, 0, -H, -1, -H, 0, H]
SIN = [0, H, 1, H, 0, -H, -1, -H]
BASE = np.array([(0.0,0.0), (1.0,0.0), (0.0,1.0)])
AXES = np.array([(math.cos(k*math.pi/4), math.sin(k*math.pi/4)) for k in range(8)])
CORNERS = [(-1.0,0.0), (1.0,0.0), (0.0,1.0)]

def sieve(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)

def lpf_d(P):
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0/P.astype(float)))[:-1]])
    return np.exp(logfac)/P

def state_arrays(P, s):
    """Return (i_s, sigma_s) where sigma_s[j]=d_{P[j]}/R_<s (defined for j>=i_s);
    the state region is area 1 so sigma_s over primes>=s telescopes to 1."""
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0/P.astype(float)))[:-1]])
    R = np.exp(logfac); d = R/P
    i_s = int(np.searchsorted(P, s))
    return i_s, d/R[i_s]

def head_count(P, s, design, p0):
    """Number of pieces in the state-s head (crossover); self-contained mirror of
    triangle_certify.tail so the builder needs no target-n."""
    logfac = np.concatenate([[0.0], np.cumsum(np.log1p(-1.0/P.astype(float)))[:-1]])
    R = np.exp(logfac); sq = np.sqrt(R/P)
    Tc = np.concatenate([[0.0], np.cumsum(sq)])[:-1]
    i_s = int(np.searchsorted(P, s)); Rs = R[i_s]; sRs = math.sqrt(Rs)
    j = np.arange(i_s, len(P)); KAP = 3.0/math.sqrt(2.0)
    if design == "uniform8":
        num = 2.0*sRs + 2.0*KAP*(Tc[j]-Tc[i_s])
    else:
        sqf = np.where(P <= p0, sq, 0.0); sqa = np.where(P > p0, sq, 0.0)
        Cf = np.concatenate([[0.0], np.cumsum(sqf)])[:-1]
        Ca = np.concatenate([[0.0], np.cumsum(sqa)])[:-1]
        num = 3.0*sRs + 2.0*(2.0*(Cf[j]-Cf[i_s]) + 1.5*(Ca[j]-Ca[i_s]))
    crit = num/np.sqrt(P[j]*R[j]) + (j-i_s)/P[j]
    bad = np.where(crit >= 1.0)[0]
    return int(bad[-1]+1) if len(bad) else 0

def templates():
    out = []
    for k in range(8):
        c, s = COS[k], SIN[k]
        v = BASE @ np.array([[c, -s], [s, c]]).T
        for a in range(3):
            out.append((k, a, v - v[a]))
    return out

TMPL = templates()
# per-template, per-vertex exact-zero flags for the three container constraints
def _flags():
    fl = []
    for (k, a, U) in TMPL:
        f = []
        for v in range(3):
            tx, ty = U[v]
            f.append((ty == 0.0, (tx + ty) == 0.0, (ty - tx) == 0.0))
        fl.append(f)
    return fl
TFLAGS = _flags()

class Builder:
    def __init__(self, eps=1e-9, cell=0.02, seed=7):
        self.eps = eps; self.cell = cell
        self.rng = np.random.default_rng(seed)
        self.placed = []; self.verts = []; self.grid = {}
        self.anchors = list(CORNERS); self.fresh = []
        self.order = None; self.fails = {}; self.jitter = 0.0

    def in_container(self, V, ti, exactL=False):
        F = TFLAGS[ti]
        for v in range(3):
            e_y, e_l, e_r = F[v]
            x, y = V[v]
            if y < (0.0 if (e_y or exactL) else self.eps): return False
            if 1.0 - x - y < (0.0 if (e_l or exactL) else self.eps): return False
            if 1.0 - y + x < (0.0 if (e_r or exactL) else self.eps): return False
        return True

    def cells_of(self, lo, hi):
        for cx in range(int(math.floor(lo[0]/self.cell)), int(math.floor(hi[0]/self.cell))+1):
            for cy in range(int(math.floor(lo[1]/self.cell)), int(math.floor(hi[1]/self.cell))+1):
                yield (cx, cy)

    def neighbors(self, lo, hi):
        out = set()
        for c in self.cells_of(lo, hi): out.update(self.grid.get(c, ()))
        return out

    def min_gap(self, V, nb):
        """Worst neighbour SAT gap: >=0 means disjoint interiors (touch allowed)."""
        if not nb: return math.inf
        N = np.stack([self.verts[i] for i in nb])
        pc = V @ AXES.T; pn = N @ AXES.T
        cmax, cmin = pc.max(0), pc.min(0)
        nmax, nmin = pn.max(1), pn.min(1)
        g = np.maximum(nmin - cmax, cmin - nmax)   # per neighbour per axis
        return float(g.max(axis=1).min())

    def contact(self, V, nb):
        """Deepest-nestle score: shared-boundary length with container + neighbours."""
        sc = 0.0; tol = 6*ETA
        E = [(V[0],V[1]), (V[1],V[2]), (V[2],V[0])]
        for (p, q) in E:
            l = math.hypot(q[0]-p[0], q[1]-p[1])
            if p[1] < tol and q[1] < tol: sc += l
            if abs(1-p[0]-p[1]) < tol and abs(1-q[0]-q[1]) < tol: sc += l
            if abs(1-p[1]+p[0]) < tol and abs(1-q[1]+q[0]) < tol: sc += l
        if nb:
            for i in nb:
                W = self.verts[i]
                for (p, q) in E:
                    d = (q[0]-p[0], q[1]-p[1]); l = math.hypot(*d)
                    u = (d[0]/l, d[1]/l); nrm = (-u[1], u[0])
                    off = p[0]*nrm[0] + p[1]*nrm[1]
                    offs = W @ np.array(nrm)
                    if np.all(np.abs(offs - off) < tol):
                        t0, t1 = sorted((0.0, l))
                        tw = np.sort((W - p) @ np.array(u))
                        ov = min(t1, tw[-1]) - max(t0, tw[0])
                        if ov > 0: sc += ov
        return sc

    def try_piece(self, L, tmpl_ids, scan_limit=600, fail_cap=48, want=10):
        if self.order is None or len(self.order) < len(self.anchors):
            pr = np.array([a[1] + 0.3*abs(a[0]) for a in self.anchors])   # bottom-up fill
            if self.jitter: pr = pr + self.jitter*self.rng.random(len(pr))
            self.order = list(np.argsort(pr))
        best = None; found = 0; seen = 0
        for idx in list(reversed(self.fresh[-48:])) + self.order:
            if self.fails.get(idx, 0) > fail_cap: continue
            ax, ay = self.anchors[idx]; hit = False
            for t in tmpl_ids:
                V = np.array([ax, ay]) + L*TMPL[t][2]
                if not self.in_container(V, t): continue
                lo, hi = V.min(0), V.max(0)
                nb = self.neighbors(lo, hi)
                if self.min_gap(V, nb) < -1e-12: continue      # touch allowed
                got = self.nudge(ax, ay, t, L)
                if got is None: continue
                sc = self.contact(got[2], self.neighbors(got[2].min(0), got[2].max(0)))
                if best is None or sc > best[0]:
                    best = (sc, t) + got
                hit = True; found += 1
            if not hit:
                self.fails[idx] = self.fails.get(idx, 0) + 1
            seen += 1
            if found >= want or seen >= scan_limit: break
        if best is None: return None
        _, t, nax, nay, V = best
        k, a, _ = TMPL[t]
        return (k, a, nax, nay, V)

    def nudge(self, ax, ay, t, L, tries=(0.0, 1.0, 2.5)):
        """Shift the anchor by ~ETA along one of the 8 axis directions (or none)
        until every neighbour gap is >= ETA, preserving containment. Snapped."""
        for m in (-1,) + tuple(range(8)):
            for mult in tries:
                if m == -1 and mult > 0: continue
                dx = 0.0 if m == -1 else ETA*mult*AXES[m][0]
                dy = 0.0 if m == -1 else ETA*mult*AXES[m][1]
                nax = round((ax+dx)/SNAP)*SNAP; nay = round((ay+dy)/SNAP)*SNAP
                V = np.array([nax, nay]) + L*TMPL[t][2]
                if not self.in_container(V, t): continue
                if self.min_gap(V, self.neighbors(V.min(0), V.max(0))) >= ETA:
                    return (nax, nay, V)
        return None

    def commit(self, p, k, a, ax, ay, V):
        i = len(self.placed)
        self.placed.append((int(p), k, a, float(ax), float(ay)))
        self.verts.append(V)
        lo, hi = V.min(0), V.max(0)
        for c in self.cells_of(lo, hi): self.grid.setdefault(c, []).append(i)
        base = len(self.anchors)
        for (vx, vy) in V:
            self.anchors.append((round(vx/SNAP)*SNAP, round(vy/SNAP)*SNAP))
        self.fresh += [base, base+1, base+2]
        if i % 200 == 0: self.order = None

    def pop_last(self, m):
        for _ in range(min(m, len(self.placed) - 1)):     # never pop the piece-2 half
            self.placed.pop(); V = self.verts.pop()
            i = len(self.placed)
            lo, hi = V.min(0), V.max(0)
            for c in self.cells_of(lo, hi):
                if i in self.grid.get(c, ()): self.grid[c].remove(i)
            del self.anchors[-3:]
        self.fresh = [f for f in self.fresh if f < len(self.anchors)]
        self.order = None; self.fails = {}

def tmpl_ids_for(design, p, p0, parity):
    if design == "uniform8" or p <= p0: return list(range(24))
    k = 0 if parity == 0 else 4
    return [t for t in range(24) if TMPL[t][0] == k]

def save(path, B, design, p0, state=2, note=""):
    doc = {"format": "triangle-head", "frame": "hypotenuse-on-x-axis",
           "design": design, "p0": p0, "state": state, "n": len(B.placed),
           "q_last": B.placed[-1][0] if B.placed else None,
           "note": note or "BUILT (verify with triangle_certify.py certify)",
           "pieces": [list(r) for r in B.placed]}
    json.dump(doc, open(path, "w"))
    print(f"  wrote {path}  (n={doc['n']}, q_last={doc['q_last']})")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-n", type=int, default=0)   # 0 -> auto from crossover
    ap.add_argument("--state", type=int, default=2)
    ap.add_argument("--design", default="uniform8", choices=["uniform8", "disciplined"])
    ap.add_argument("--p0", type=int, default=1000)
    ap.add_argument("--out", default="triangle_head.json")
    ap.add_argument("--checkpoint", default=""); ap.add_argument("--ckpt-every", type=int, default=500)
    ap.add_argument("--resume", default="")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--time-cap", type=float, default=0.0)
    a = ap.parse_args()

    s = a.state
    P = sieve(max(1_200_000, 20*max(a.target_n, 1)))
    # auto target from the crossover if not given
    tgt = a.target_n if a.target_n > 0 else head_count(P, s, a.design, a.p0)
    if tgt == 0:
        print(f"state {s} is terminal (empty head); nothing to build."); return
    P = sieve(max(1_200_000, 20*tgt))
    i_s, sig = state_arrays(P, s)                       # sig[j] = sigma_s(P[j]), j>=i_s
    Pq = P[i_s:]; Lq = np.sqrt(2.0*sig[i_s:])           # state-local prime & leg arrays
    print(f"building state s={s}, design={a.design}, target n={tgt} "
          f"(primes {int(Pq[0])}..{int(Pq[tgt-1])})")
    B = Builder(seed=a.seed); start = 0
    if a.resume:
        doc = json.load(open(a.resume))
        for (p, k, aa, ax, ay) in doc["pieces"]:
            li = int(np.searchsorted(Pq, p))
            t = [t_ for t_ in range(24) if TMPL[t_][0]==k and TMPL[t_][1]==aa][0]
            B.commit(p, k, aa, ax, ay, np.array([ax, ay]) + Lq[li]*TMPL[t][2])
        start = len(B.placed); print(f"resumed {start} pieces (q_last={doc['q_last']})")
    if start == 0 and s == 2:
        B.commit(2, 0, 0, 0.0, 0.0, np.array([0.0,0.0]) + Lq[0]*TMPL[0][2])  # exact half
        start = 1

    t0 = time.time(); backoffs = [25, 100, 400, 1200]
    i = start
    while i < tgt:
        p = int(Pq[i]); ids = tmpl_ids_for(a.design, p, a.p0, i % 2)
        got = B.try_piece(Lq[i], ids)
        if got is None:
            B.fails = {}; B.order = None
            got = B.try_piece(Lq[i], ids, scan_limit=10**9, fail_cap=10**9, want=1)
        if got is None:
            repaired = False
            for B_off in backoffs:
                print(f"  jam at piece {i} (p={p}); backing off {B_off}")
                B.pop_last(B_off); B.jitter = 0.02*float(B.rng.random())
                j = max(len(B.placed), 1 if s == 2 else 0); ok = True
                while j <= i:
                    pj = int(Pq[j]); idj = tmpl_ids_for(a.design, pj, a.p0, j % 2)
                    g = B.try_piece(Lq[j], idj, scan_limit=10**9, fail_cap=10**9)
                    if g is None: ok = False; break
                    B.commit(pj, *g); j += 1
                if ok: repaired = True; B.jitter = 0.0; break
            if not repaired:
                print(f"UNREPAIRED JAM at piece {i} (p={p}); saving partial head")
                save(a.out, B, a.design, a.p0, s, note=f"PARTIAL: jammed at {p}"); return
            i = len(B.placed); continue
        B.commit(p, *got); i += 1
        if i % 250 == 0:
            print(f"  {i}/{tgt} pieces (p={p}), covered {float(sig[i_s:i_s+i].sum()):.5f}, "
                  f"{time.time()-t0:.0f}s")
        if a.checkpoint and i % a.ckpt_every == 0:
            save(a.checkpoint, B, a.design, a.p0, s, note="CHECKPOINT")
        if a.time_cap and time.time()-t0 > a.time_cap:
            print("time cap reached"); break
    save(a.out, B, a.design, a.p0, s)
    print(f"done: state {s}, {len(B.placed)} pieces in {time.time()-t0:.0f}s, "
          f"covered {float(sig[i_s:i_s+len(B.placed)].sum()):.6f}")

if __name__ == "__main__":
    main()
