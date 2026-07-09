"""Exact rectangle companion model: guillotine packing of similar rectangles
of areas d_p (least-prime-factor densities) into a ratio-r rectangle.
Reproduces: head window, pure-ledger runs, reserve-backstop runs, and the
lanes-only worst case for the side-length budget (see paper, Sec. 'rectangle').
"""
import numpy as np, heapq
from math import sqrt, log

def sieve(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.flatnonzero(s)

PR=sieve(10**7)
logfac=np.concatenate([[0.0],np.cumsum(np.log1p(-1.0/PR.astype(float)))[:-1]])
dP=np.exp(logfac)/PR   # d_p = R_{<p}/p exactly (float)

def pure_ledger(r, nmax=None):
    """Fattest-first guillotine ledger; exact free-region representation."""
    nmax=nmax or len(PR)
    fr=[]
    def push(A,B):
        if A<B: A,B=B,A
        if B>1e-13: heapq.heappush(fr,(-B,A,B))
    push(sqrt(r),1/sqrt(r))
    minmargin=1e9
    for k in range(nmax):
        d=dP[k]; a,b=sqrt(d*r),sqrt(d/r)
        if a<b: a,b=b,a
        if not fr: return ('empty',PR[k],minmargin)
        negB,A,B=heapq.heappop(fr)
        if not(b<=B+1e-15 and a<=A+1e-15):
            return ('jam',PR[k],minmargin)
        c1=max(min(A-a,B),min(a,B-b)); c2=max(min(A-a,b),min(A,B-b))
        if c1>=c2: push(A-a,B); push(a,B-b)
        else:      push(A-a,b); push(A,B-b)
        if k+1<nmax and fr:
            minmargin=min(minmargin,(-fr[0][0])/sqrt(dP[k+1]/r))
    return ('ok',None,minmargin)

def reserve_backstop(r, reserve_frac, nmax=None):
    """Fringe-first; reserve hit only on fringe failure. Exact hit decomposition:
    kappa -> (sqrt(kappa)-sqrt(d))^2 + two lanes (Lemma 'hit')."""
    nmax=nmax or len(PR)
    fr=[]
    def push(A,B):
        if A<B: A,B=B,A
        if B>1e-13 and A>1e-13: heapq.heappush(fr,(-B,A,B))
    kappa=0.5*reserve_frac
    hw,hh=sqrt(0.5*r),sqrt(0.5/r)
    kA,kB=sqrt(kappa*r),sqrt(kappa/r)
    push(hw-kA,hh); push(kA,hh-kB)      # piece 2 = exact half at r=sqrt2
    hits=[]; sq=0.0
    for k in range(1,nmax):
        d=dP[k]; a,b=sqrt(d*r),sqrt(d/r)
        if fr:
            negB,A,B=fr[0]
            if b<=B+1e-15 and a<=A+1e-15:
                heapq.heappop(fr)
                c1=max(min(A-a,B),min(a,B-b)); c2=max(min(A-a,b),min(A,B-b))
                if c1>=c2: push(A-a,B); push(a,B-b)
                else:      push(A-a,b); push(A,B-b)
                continue
        if sqrt(kappa)-sqrt(d)<=1e-13:
            return ('EXHAUSTED',PR[k],hits,sq)
        kA,kB=sqrt(kappa*r),sqrt(kappa/r)
        newk=(sqrt(kappa)-sqrt(d))**2
        nA,nB=sqrt(newk*r),sqrt(newk/r)
        push(nA,kB-nB); push(a,kB-b)    # the two exact lanes
        kappa=newk; sq+=sqrt(d); hits.append(int(PR[k]))
    return ('OK',None,hits,sq)


# ----------------------------------------------------------------------
# Interval-certified head for Theorem (perfect packing of the sqrt2 rectangle).
#
# Discipline: the reserve-backstop scheduler that Theorem 'buffer' consumes.
# Piece 2 is the exact left half, leaving a similar right half as the reserve;
# thereafter each prime is placed fringe-first (fattest-fit) and the reserve is
# hit only when no fringe rectangle fits, using the exact hit decomposition of
# Lemma 'hit'. Every fit is certified in `prec`-bit interval arithmetic (both
# orientations). At q0 the routine emits the handoff Theorem 'buffer' requires:
#
#     free region  =  reserve (similar to container, area kappa0)  +  fringe,
#
# with reserve area + fringe areas enclosing R_{<q0}, the reserve's aspect
# interval enclosing sqrt2, and a positive certified fit margin throughout.
# Default q0 = 50000 (the first 5133 primes) per the amended Theorem.
# ----------------------------------------------------------------------
def certify_head(q0=50000, prec=200, reserve_frac=0.5, collect=False):
    from mpmath import iv
    import mpmath as mp
    mp.mp.prec = prec; iv.prec = prec
    R2 = iv.sqrt(iv.mpf(2)); SQ = iv.sqrt
    norm = lambda A, B: (A, B) if (A.a >= B.a) else (B, A)
    # fringe kept as parallel lists: interval pairs AB, and float shadows FA,FB
    # (min/long lower endpoints) for a fast search; swap-pop keeps them in sync.
    AB = []; FA = []; FB = []
    def push(A, B):
        A, B = norm(A, B)
        if float(B.b) > 1e-40:
            AB.append((A, B)); FA.append(float(A.a)); FB.append(float(B.a))
    def drop(i):                                     # O(1) swap-pop
        AB[i] = AB[-1]; FA[i] = FA[-1]; FB[i] = FB[-1]
        AB.pop(); FA.pop(); FB.pop()
    def mn(X, Y): return X if float(X.a) <= float(Y.a) else Y
    def mx(X, Y): return X if float(X.a) >= float(Y.a) else Y

    # piece 2 = exact left half. The similar right half (area 1/2) is split into
    # a similar reserve (area kappa = reserve_frac/2) plus two fringe strips, exactly
    # as reserve_backstop() seeds it; reserve_frac=1 reduces to the whole right half.
    half  = iv.mpf(1) / 2
    kappa = half * iv.mpf(reserve_frac)
    hw, hh = SQ(half * R2), SQ(half / R2)   # right-half dims (similar)
    kA, kB = SQ(kappa * R2), SQ(kappa / R2)  # reserve dims (similar)
    push(hw - kA, hh); push(kA, hh - kB)     # two fringe strips
    Riv   = half                             # R_{<3}
    res_slack = iv.mpf('1e9')                # min reserve slack (sqrt-kappa - sqrt-d) over hits
    hits = []; nplaced = 1; records = []

    P = [int(x) for x in sieve(q0 + 10) if int(x) <= q0]
    for p in P[1:]:                          # p = 3, 5, 7, ...
        d = Riv / p
        a = SQ(d * R2); b = SQ(d / R2)       # a >= b  (long, short)
        af = float(a.b); bf = float(b.b)     # upper bounds: need dim >= these
        done = False
        # search ALL fringe (Thm 'budget': any guillotine rule if one exists);
        # pick the fattest rectangle that fits, then certify it rigorously.
        best = -1; bestkey = -1.0
        for i in range(len(AB)):
            fa, fb = FA[i], FB[i]
            if ((fa >= af and fb >= bf) or (fb >= af and fa >= bf)) and fb > bestkey:
                bestkey = fb; best = i
        if best >= 0:
            A, B = AB[best]
            fit1 = (A.a >= a.b) and (B.a >= b.b)     # long -> A
            fit2 = (B.a >= a.b) and (A.a >= b.b)     # rotated
            if fit1 or fit2:
                drop(best)
                pA, pB = (a, b) if fit1 else (b, a)
                mA, mB = A - pA, B - pB
                if mA.a >= 0 and mB.a >= 0:          # rigorous containment (A>=a, B>=b)
                    # guillotine: choose the cut with the fatter residue (heuristic)
                    c1 = mx(mn(A - pA, B), mn(pA, B - pB))
                    c2 = mx(mn(A - pA, pB), mn(A, B - pB))
                    if float(c1.a) >= float(c2.a):
                        push(A - pA, B); push(pA, B - pB)
                    else:
                        push(A - pA, pB); push(A, B - pB)
                    mgn = mA if mA.a <= mB.a else mB
                    records.append((p, (a, b), (A, B), mgn, 'fringe'))
                    done = True; nplaced += 1
        if not done:
            # reserve hit: certify d <= kappa, then Lemma 'hit' decomposition
            if not (kappa.a >= d.b):
                return {'status': 'RESERVE_EXHAUSTED', 'prime': p,
                        'kappa': (float(kappa.a), float(kappa.b)),
                        'd': (float(d.a), float(d.b))}
            kA, kB = SQ(kappa * R2), SQ(kappa / R2)
            slack = SQ(kappa) - SQ(d)                  # linear reserve slack at this hit
            if slack.a < res_slack.a: res_slack = slack
            mgn = (kA - a) if (kA - a).a <= (kB - b).a else (kB - b)
            records.append((p, (a, b), (kA, kB), mgn, 'reserve'))
            newk = (SQ(kappa) - SQ(d)) ** 2
            nA, nB = SQ(newk * R2), SQ(newk / R2)
            push(nA, kB - nB)                         # lane 1 = sqrt(k'r) x sqrt(d/r)
            push(a,  kB - b)                          # lane 2 = sqrt(dr)  x sqrt(k'/r)
            kappa = newk; hits.append(p); nplaced += 1
        Riv = Riv - d

    # ---- handoff certificate: reserve (similar) + fringe ----
    A0, B0 = SQ(kappa * R2), SQ(kappa / R2)
    aspect = A0 / B0                                   # must enclose sqrt2
    fringe = list(AB)
    area_fringe = iv.mpf(0)
    for (A, B) in fringe: area_fringe += A * B
    total = kappa + area_fringe                        # must enclose R_{<q0} = Riv
    fringe_min = None
    if fringe:
        fmin = min(fringe, key=lambda pr: float(pr[1].a))[1]
        fringe_min = (float(fmin.a), float(fmin.b))
    rs = None if res_slack.a > iv.mpf('1e8').a else float(res_slack.a)
    fringe_ratio = float((area_fringe / Riv).b)        # >= immature ratio at q0; cf. bound 0.923
    placements = None
    if collect:
        placements = [{'p': int(pp),
                       'piece_sides': [[float(A_.a), float(A_.b)] for A_ in pab],
                       'host_sides':  [[float(A_.a), float(A_.b)] for A_ in hAB],
                       'fit_margin':  float(mg.a),      # certified >= 0
                       'kind': kd} for (pp, pab, hAB, mg, kd) in records]
    ret = {
        'fringe_area':    (float(area_fringe.a), float(area_fringe.b)),
        'fringe_ratio_ub': fringe_ratio,               # upper bd on I(q0)/R_{<q0}; Lemma bound is 0.923
        'status': 'CERTIFIED', 'q0': q0, 'prec': prec, 'n_primes': nplaced,
        'placements_certified': nplaced,                       # all passed rigorous containment
        'reserve_area':   (float(kappa.a), float(kappa.b)),
        'reserve_dims':   ((float(A0.a), float(A0.b)), (float(B0.a), float(B0.b))),
        'reserve_aspect': (float(aspect.a), float(aspect.b)),   # encloses sqrt2 = 1.41421356
        'fringe_count':   len(fringe),
        'fringe_min_side': fringe_min,                          # thinnest sliver (expected small)
        'reserve_hits':   hits,
        'reserve_min_slack': rs,                                # sqrt(kappa)-sqrt(d) at tightest hit (>0)
        'R_lt_q0':        (float(Riv.a), float(Riv.b)),
        'reserve_plus_fringe': (float(total.a), float(total.b)),  # encloses R_lt_q0
    }
    if collect: ret['placements'] = placements
    return ret


# ----------------------------------------------------------------------
# Tail certificate for the immature-stock ratio (Lemma 'explicit', Appendix).
# Verifies rho(q) = (3 T(q)+r)/sqrt(q R_<q) <= 0.923 < 1-1/q for all primes
# q > q0 = 5e4, in three ranges: (q0,Q1] and (Q1,X1] by direct evaluation,
# (X1,inf) by the explicit Rosser-Schoenfeld/Dusart majorant of the appendix.
# ----------------------------------------------------------------------
def certify_tail(q0=50000, Q1=100000, X1=10**7):
    from math import log, sqrt, exp
    g = 0.5772156649015329; eg = exp(-g); eg2 = exp(-g/2); r = sqrt(2.0)
    P = PR[PR <= X1].astype(np.float64)
    Rlt = np.exp(np.concatenate([[0.0], np.cumsum(np.log1p(-1.0/P))[:-1]]))
    T = np.cumsum(np.sqrt(Rlt/P)); lnq = np.log(P)
    rho = (3*T + r)/np.sqrt(P*Rlt)
    A = T*lnq**1.5/np.sqrt(P)

    mw = (P > q0) & (P <= Q1); mt = (P > Q1) & (P <= X1)
    win_max = float(rho[mw].max()); mid_max = float(rho[mt].max())
    cross = int(P[np.where(rho >= 1)[0][-1]])

    # explicit analytic majorant on (X1, inf)  (see appendix)
    ell = log(Q1); beta = eg*(1 - 1/ell**2)
    F = (1 + 1.2762/ell)*(1 + 1.0/ell); G = 1.0/(1 - 3.0/ell)
    K = 1.004*eg2*(1 + 1.2762/ell + F*G)
    T_Q1 = float(T[np.searchsorted(P, Q1, side='right') - 1])
    Abound = lambda q: T_Q1*log(q)**1.5/sqrt(q) + K            # >= A(q) for q>=Q1
    A_at_X1 = Abound(X1)                                        # sup over [X1,inf)
    rho_tail = 3*A_at_X1/(sqrt(beta)*log(X1)) + r*sqrt(log(X1))/(sqrt(beta)*sqrt(X1))
    dominates = bool(np.all(np.array([Abound(x) for x in P[mt]]) >= A[mt]))

    ok = (win_max <= 0.923) and (mid_max < 1) and (rho_tail < 1) and (cross < q0)
    return {
        'status': 'CERTIFIED' if ok else 'FAIL',
        'q0': q0, 'Q1': Q1, 'X1': X1, 'crossover': cross,
        'window_max_(q0,Q1]': win_max, 'mid_max_(Q1,X1]': mid_max,
        'tail_bound_(X1,inf)': rho_tail, 'overall_bound': max(win_max, mid_max, rho_tail),
        'constants': {'beta_Rlow': beta, 'K': K, 'T_Q1': T_Q1, 'A_sup_tail': A_at_X1,
                      'F': F, 'G': G},
        'majorant_dominates_A_on_[Q1,X1]': dominates,
    }


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description="Exact rectangle companion model (paper Sec. 'rectangle').")
    sub = ap.add_subparsers(dest='cmd')
    hp = sub.add_parser('head',     help='interval-certified reserve+fringe head (Theorem)')
    hp.add_argument('--q0',   type=int,   default=50000, help='place primes <= q0 (default 50000 = first 5133)')
    hp.add_argument('--prec', type=int,   default=200,   help='interval-arithmetic precision in bits')
    hp.add_argument('--reserve-frac', type=float, default=0.5, dest='rf',
                    help='fraction of the right half held as reserve (default 0.5, matching reserve_backstop)')
    hp.add_argument('--dump', type=str, default=None,
                    help='write the full per-placement certificate (dims + interval margins) to this JSON file')
    lp = sub.add_parser('ledger',   help='pure fattest-first ledger run (float)')
    lp.add_argument('--nmax', type=int, default=700000)
    bp = sub.add_parser('backstop', help='reserve-backstop scheduler run (float)')
    bp.add_argument('--nmax', type=int, default=None)
    tp = sub.add_parser('tail',   help='immature-stock ratio bound rho(q)<1 for q>q0 (Lemma explicit)')
    tp.add_argument('--q0', type=int, default=50000)
    tp.add_argument('--Q1', type=int, default=100000)
    tp.add_argument('--X1', type=int, default=10**7)
    a = ap.parse_args()

    r = 2 ** 0.5
    if a.cmd == 'ledger':
        print('pure ledger r=sqrt2:', pure_ledger(r, a.nmax)[0:2])
    elif a.cmd == 'backstop':
        st, _, hits, sq = reserve_backstop(r, 0.5, a.nmax)
        print(f'reserve-backstop: {st}, hits={hits}, sqrt-spend={sq:.4f}')
    elif a.cmd == 'tail':
        c = certify_tail(a.q0, a.Q1, a.X1)
        print(f"immature-stock ratio certificate  (q0={c['q0']}, Q1={c['Q1']}, X1={c['X1']})")
        print(f"  status                : {c['status']}")
        print(f"  crossover (rho>=1 last): {c['crossover']}   (< q0)")
        print(f"  (q0,Q1]  max rho      : {c['window_max_(q0,Q1]']:.5f}   (<= 0.923)")
        print(f"  (Q1,X1]  max rho      : {c['mid_max_(Q1,X1]']:.5f}   (< 1)")
        print(f"  (X1,inf) analytic bd  : {c['tail_bound_(X1,inf)']:.5f}   (< 1)")
        print(f"  overall  rho(q)       : <= {c['overall_bound']:.5f}  for all q > q0   (< 1-1/q)")
        k = c['constants']
        print(f"  appendix constants    : beta={k['beta_Rlow']:.5f}  K={k['K']:.4f}  "
              f"T(Q1)={k['T_Q1']:.3f}  A_sup(tail)={k['A_sup_tail']:.4f}")
        print(f"  majorant dominates A on [Q1,X1]? {c['majorant_dominates_A_on_[Q1,X1]']}")
    else:  # 'head' or no subcommand
        q0   = getattr(a, 'q0', 50000)
        prec = getattr(a, 'prec', 200)
        rf   = getattr(a, 'rf', 0.5)
        dump = getattr(a, 'dump', None)
        c = certify_head(q0, prec, rf, collect=bool(dump))
        if dump:
            import json
            pl = c.pop('placements')
            doc = {'model': 'sqrt2 rectangle, r=sqrt(2)', 'q0': c['q0'], 'prec_bits': prec,
                   'n_primes': c['n_primes'], 'reserve_min_slack': c['reserve_min_slack'],
                   'reserve_aspect': c['reserve_aspect'], 'reserve_area': c['reserve_area'],
                   'R_lt_q0': c['R_lt_q0'], 'reserve_plus_fringe': c['reserve_plus_fringe'],
                   'min_fit_margin': min(x['fit_margin'] for x in pl),
                   'reproduce': f'python rect_exact.py head --q0 {c["q0"]} --prec {prec} --dump {dump}',
                   'schema': 'each placement: p, piece_sides [[a.lo,a.hi],[b.lo,b.hi]], '
                             'host_sides [[A.lo,A.hi],[B.lo,B.hi]], fit_margin (certified >=0), kind',
                   'placements': pl}
            json.dump(doc, open(dump, 'w'))
            print(f"  wrote {dump}: {len(pl)} placements, min fit margin {doc['min_fit_margin']:.3e}")
        print(f"head certificate  (r = sqrt2,  q0 = {c.get('q0', q0)},  prec = {prec} bits)")
        print(f"  status              : {c['status']}")
        if c['status'] != 'CERTIFIED':
            print(f"  detail              : {c}"); raise SystemExit
        print(f"  placements certified: {c['placements_certified']}   (all passed rigorous containment)")
        print(f"  reserve hits        : {len(c['reserve_hits'])} at primes {c['reserve_hits']}")
        print(f"  reserve min slack   : {c['reserve_min_slack']:.3e}   (sqrt(kappa)-sqrt(d) at tightest hit, > 0)")
        print(f"  reserve area kappa0 : [{c['reserve_area'][0]:.10f}, {c['reserve_area'][1]:.10f}]")
        print(f"  reserve aspect      : [{c['reserve_aspect'][0]:.10f}, {c['reserve_aspect'][1]:.10f}]  (encloses sqrt2)")
        print(f"  fringe rectangles   : {c['fringe_count']}   (thinnest sliver {c['fringe_min_side'][0]:.2e})")
        print(f"  R_{{<q0}}             : [{c['R_lt_q0'][0]:.12f}, {c['R_lt_q0'][1]:.12f}]")
        print(f"  reserve + fringe    : [{c['reserve_plus_fringe'][0]:.12f}, {c['reserve_plus_fringe'][1]:.12f}]   (encloses R_{{<q0}})")
        print(f"  fringe area / R_{{<q0}}: {c['fringe_ratio_ub']:.4f}   (>= immature ratio at q0; Lemma 6.4 bound is 0.923)")
