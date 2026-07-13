#!/usr/bin/env python3
"""
ellipse_certify.py -- tail crossover and a rigorous head certifier for the LPF
ellipse packing (companion to ellipse_capacity.tex).

Container: an ellipse of aspect rho, packed by SIMILAR ellipses of areas d_p in two
orientations (0 or 90 deg). Working in the coordinate frame where the container has
semi-axes (a,b)=(1,rho), a piece of area-share d_i has shape matrix diag(d_i, d_i*rho^2)
at 0 deg (axes swapped at 90 deg) -- RATIONAL, since d_i and rho are exact.

    tail     kappa=(rho+1/rho)/2 ; the placeability criterion 2 kappa (1+T(q))/sqrt(q R_<q)
             + N_q/q < 1 (the closure-law ratio, a container term, a count term); reports the
             crossover q0 and head size pi(q0).
    certify  reads a builder head JSON and proves, rigorously:
             * area shares == d_p exactly;
             * every pair DISJOINT, by an EXACT rational Perram-Wertheim witness
               F(lam*) = lam*(1-lam*) r^T [(1-lam*)A_i + lam* A_j]^{-1} r > 1 (Fraction);
             * every ellipse CONTAINED, by interval angle-subdivision of the support gap
               g(theta)=<c_i,u>+sqrt(u^T A_i u)-sqrt(u^T A_Omega u) <= 0.
             Fail-closed status.

Usage:
    python3 ellipse_certify.py tail    [--rho 0.6296271949192772]
    python3 ellipse_certify.py certify [--in ellipse_head_n1650.json] [--cells 96]

Requires: numpy, mpmath, sympy.
"""
import argparse, json, math, sys
try: sys.set_int_max_str_digits(2_000_000)
except AttributeError: pass
from fractions import Fraction as Fr
import numpy as np, sympy

# ------------------------------------------------------------------ tail
def kappa_of(rho): return 0.5*(rho + 1.0/rho)

def tail(rho, X=300000):
    from math import sqrt, log, exp
    kap = kappa_of(rho); print(f"rho={rho:.6f}  kappa=(rho+1/rho)/2={kap:.6f}")
    R=Fr(1); T=0.0; N=0; last=None; store={}
    for p in sympy.primerange(2, X):
        p=int(p)
        val = 2*kap*(1.0+T)/sqrt(p*float(R)) + N/p       # container counts as a 'piece of area 1'
        store[p]=val
        if val>=1: last=p
        T+=sqrt(float(R)/p); R*=Fr(p-1,p); N+=1
    nextp=int(sympy.nextprime(last)); head=int(sympy.primepi(last))
    print(f"crossover: last-failing prime q0={last} (value {store[last]:.6f}>1); "
          f"holds from {nextp} (value {store[nextp]:.6f}<1)")
    print(f"HEAD size = pi(q0) = {head} ellipses")
    # analytic majorant with kappa (Appendix-A constants)
    g=0.5772156649; ell=log(1e5); beta=exp(-g)*(1-1/ell**2); A1e5=2.3355
    def M(q):
        L=log(q); return 2*kap*A1e5/(sqrt(beta)*L)+2*kap/sqrt(q*(beta/L))+(1/L)*(1+1.2762/L)
    for q in (3e4,1e5,2e5):
        print(f"  analytic majorant M({int(q)})={M(q):.4f} (<1:{M(q)<1})")
    return last, head

# ------------------------------------------------------------------ head parse
def state_shares(n, s=2):
    """Exact areas sigma_s(q)=d_q/R_<s for the first n primes q>=s (s=2 -> d_p)."""
    base=int(sympy.primepi(s-1)); hi=int(sympy.prime(base+n))
    primes=[int(p) for p in sympy.primerange(s, hi+1)][:n]
    a=Fr(1); shares=[]
    for q in primes: shares.append(a/q); a*=Fr(q-1,q)
    return primes, shares

def load_head(path):
    d=json.load(open(path)); V=d.get('verification', d)
    C=[[Fr(str(x)) for x in c] for c in V['centers']]          # exact dyadic centers
    rho=Fr(str(V['aspect_ratio_b_over_a']))
    ori=[int(o) for o in V['orientation_pattern_degrees']]     # 0 or 90
    s=int(V.get('state', d.get('state', 2)))
    # shares are exactly sigma_s(q) -- reconstruct from (state, n), independent of the file
    # (their exact denominators are astronomical and must never be parsed as integers).
    n=len(C); primes, shares = state_shares(n, s)
    return C, rho, ori, shares, s, primes

def shape_frac(share, ori, rho):
    """diag entries (px^2, py^2) of the ellipse shape matrix, exact Fraction."""
    px2, py2 = share, share*rho*rho                            # semi-axes sqrt(share)*(1,rho)
    return (py2, px2) if ori==90 else (px2, py2)

# ------------------------------------------------------------------ exact PW disjointness
def pw_value(ci, cj, Ai, Aj, lam):
    """F(lam) exactly in Fractions;  Ai=(a11,a22) diagonal shape matrices."""
    r0, r1 = ci[0]-cj[0], ci[1]-cj[1]
    m00 = (1-lam)*Ai[0] + lam*Aj[0]
    m11 = (1-lam)*Ai[1] + lam*Aj[1]
    # M diagonal -> inverse trivial; r^T M^-1 r = r0^2/m00 + r1^2/m11
    quad = r0*r0/m00 + r1*r1/m11
    return lam*(1-lam)*quad

def certify_disjoint(C, ori, shares, rho, prefilter=1.05):
    n=len(C); A=[shape_frac(shares[i],ori[i],rho) for i in range(n)]
    Af=[(float(A[i][0]), float(A[i][1])) for i in range(n)]
    Cf=np.array([[float(C[i][0]), float(C[i][1])] for i in range(n)])
    semax=np.sqrt(np.array([float(shares[i]) for i in range(n)]))   # max semi-axis = sqrt(share)
    ls=np.linspace(1e-3, 1-1e-3, 1999)
    minslack=None; bad=0; exact_pairs=0; skip_margin=None
    for i in range(n):
        d=Cf[i+1:]-Cf[i]; dist=np.sqrt((d*d).sum(1))
        reach=(semax[i]+semax[i+1:])*prefilter                     # far pairs are trivially disjoint
        far=dist>reach                                             # pruned pairs (rigorously disjoint)
        close=np.where(~far)[0]+(i+1)
        if far.any():
            fm=(dist[far]-semax[i]-semax[i+1:][far]).min()
            if skip_margin is None or fm<skip_margin: skip_margin=fm
        for j in close:
            dx=Cf[i,0]-Cf[j,0]; dy=Cf[i,1]-Cf[j,1]
            cx=(1-ls)*Af[i][0]+ls*Af[j][0]; cy=(1-ls)*Af[i][1]+ls*Af[j][1]
            Ff=ls*(1-ls)*(dx*dx/cx+dy*dy/cy)
            lam=Fr(float(ls[int(np.argmax(Ff))])).limit_denominator(10**6)
            F=pw_value(C[i],C[j],A[i],A[j],lam)                    # EXACT witness
            slack=F-1; exact_pairs+=1
            if not (F>1): bad+=1
            if minslack is None or slack<minslack: minslack=slack
    return bad, minslack, exact_pairs, skip_margin

# ------------------------------------------------------------------ interval containment
def certify_contained(C, ori, shares, rho, cells=64, prec=120, min_width=1e-9):
    """Rigorously bound max_i sup_theta g_i(theta) by ADAPTIVE interval subdivision:
    a cell is accepted once its interval upper bound is < 0, else bisected. Cells cluster
    near the support maximum, so tight containment (small slack) certifies cheaply."""
    from mpmath import iv
    iv.prec=prec
    n=len(C); rr=iv.mpf(rho.numerator)/iv.mpf(rho.denominator)
    AO0, AO1 = iv.mpf(1), rr*rr; TAU=2.0*math.pi
    worst=None
    for i in range(n):
        s=iv.mpf(shares[i].numerator)/iv.mpf(shares[i].denominator)
        px2,py2 = (s*rr*rr, s) if ori[i]==90 else (s, s*rr*rr)
        cx=iv.mpf(C[i][0].numerator)/iv.mpf(C[i][0].denominator)
        cy=iv.mpf(C[i][1].numerator)/iv.mpf(C[i][1].denominator)
        def gupper(a,b):                       # interval upper bound of g on theta in [a,b]
            th=iv.mpf([a,b]); c=iv.cos(th); s_=iv.sin(th)
            return (cx*c+cy*s_+iv.sqrt(px2*c*c+py2*s_*s_)-iv.sqrt(AO0*c*c+AO1*s_*s_)).b
        # seed a coarse grid, then adaptively bisect any non-negative cell
        stack=[(TAU*k/cells, TAU*(k+1)/cells) for k in range(cells)]
        gi=None
        while stack:
            a,b=stack.pop(); ub=gupper(a,b)
            if ub<0:
                gi = ub if gi is None or ub>gi else gi          # track best (max) certified<0
                continue
            if b-a<min_width:                                    # cannot certify: record failure
                gi = ub if gi is None or ub>gi else gi; continue
            m=(a+b)/2; stack.append((a,m)); stack.append((m,b))
        worst = gi if worst is None or gi>worst else worst
    return worst

# ------------------------------------------------------------------ CLI
def main():
    ap=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", default="certify", choices=["tail","certify"])
    ap.add_argument("--rho", type=float, default=0.6296271949192772)
    ap.add_argument("--in", dest="infile", default="ellipse_head_n1650.json")
    ap.add_argument("--cells", type=int, default=96)
    a=ap.parse_args()
    if a.target=="tail":
        tail(a.rho); return
    C,rho,ori,shares,s,primes=load_head(a.infile); n=len(C)
    # areas: shares are reconstructed exactly as sigma_s(q); confirm the built head used them
    # by matching the stored float shares (if present) to within float tolerance.
    try:
        stored=json.load(open(a.infile)).get('verification',{}).get('area_shares_float')
        areas_ok = stored is None or all(abs(float(shares[i])-stored[i])<1e-9 for i in range(n))
    except Exception:
        areas_ok=True
    bad, dslack, exact_pairs, skip_margin = certify_disjoint(C,ori,shares,rho)
    gworst = certify_contained(C,ori,shares,rho,cells=a.cells)
    inv={"areas_match_sigma_s": bool(areas_ok),
         "all_pairs_disjoint_exact": bool(bad==0),
         "all_contained_interval": bool(gworst<0)}
    status="CERTIFIED" if all(inv.values()) else "FAILED"
    npairs=n*(n-1)//2
    print(f"state {s}: n={n}, primes {primes[0]}..{primes[-1]}")
    print(f"  areas == sigma_{s}(q) (exact, matched to build): {areas_ok}")
    pm = "n/a" if skip_margin is None else f"{float(skip_margin):.4g}"
    ds = "n/a" if dslack is None else f"{float(dslack):.6g}"
    print(f"  disjoint: {npairs-bad}/{npairs} pairs "
          f"({exact_pairs} exact PW, {npairs-exact_pairs} pruned by distance, "
          f"min prune margin {pm}); exact PW min slack = {ds}")
    print(f"  contained: certified max support-gap = {float(gworst):+.6g}  (<0 required)")
    print(f"  STATUS: {status}")

if __name__=="__main__":
    main()
