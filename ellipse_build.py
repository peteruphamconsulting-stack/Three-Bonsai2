#!/usr/bin/env python3
"""
ellipse_build.py -- extend a certified LPF ellipse head toward the 1650-piece tail
crossover, fast enough to finish locally.

Same geometry/margins as ellipse_binary_head_builder_v2.py (centers_normalized rows (x,u),
actual centre (x, b*u); container support sqrt(cos^2+b^2 sin^2); Perram-Wertheim pair margin
max_lambda F - 1; overall insertion margin = min(containment, min_j pair)). Two changes make it
finish locally: a SPATIAL HASH plus a short big-piece list, so each candidate is tested against
O(1) neighbours instead of all n placed pieces; and a BEAM search (like your v2 builder) that
keeps several configurations and expands each by its best few *separated* placements -- greedy
largest-margin alone paints into a corner (a locally optimal placement that globally blocks the
next piece), and no amount of reseeding escapes it, whereas the beam keeps the non-greedy
alternative that leaves room. Resume from your n=16/n=21 head; only the tail is new work.

    build   resume a head/checkpoint JSON and place primes up to --target-n by beam search:
            at each step every beam state proposes its top --per-parent separated placements
            (vectorized containment + big-piece PW prescreen, then exact small-piece margin,
            then Nelder-Mead refinement, gated on the exact margin the certifier will see);
            the best --beam children are kept. Checkpoints every --checkpoint-every pieces.
            On a stall it reshuffles the Sobol field; if the whole beam is exhausted it saves a
            checkpoint you can resume with a different --seed or a larger --beam.

Output JSON matches the builder (aspect_ratio_b_over_a, orientation_pattern_degrees, centers,
centers_normalized, area_shares) so `python3 ellipse_certify.py certify --in head.json` verifies it.

Usage:
    python3 ellipse_build.py --resume ellipse_binary_checkpoint_v2_n32.json --target-n 1650 \
        --out ellipse_head_n1650.json --checkpoint ckpt.json --beam 6 --sobol 3072 --seed 7
    # if it stalls, resume from the checkpoint with a different seed / wider beam:
    python3 ellipse_build.py --resume ckpt.json --target-n 1650 --out ellipse_head_n1650.json \
        --checkpoint ckpt.json --beam 10 --sobol 4096 --seed 11
Requires: numpy, scipy, sympy.
"""
import argparse, json, math, time, os, sys
try: sys.set_int_max_str_digits(2_000_000)
except AttributeError: pass
from fractions import Fraction as Fr
import numpy as np, sympy
from scipy.stats import qmc
from scipy.optimize import minimize

# ---------------------------------------------------------------- LPF data
def lpf_state(n, s=2):
    """First n children of a state-s cell: primes q>=s with relative areas
    sigma_s(q) = (1/q) prod_{s<=r<q}(1-1/r) = d_q/R_<s.  s=2 is the container family d_p."""
    base=int(sympy.primepi(s-1)); hi=int(sympy.prime(base+n))
    primes=[int(p) for p in sympy.primerange(s, hi+1)][:n]
    a=Fr(1); shares=[]
    for q in primes: shares.append(a/q); a*=Fr(q-1,q)       # a = prod_{s<=r<q}(1-1/r)
    scales=np.sqrt(np.array([float(x) for x in shares]))
    return primes, shares, scales

def lpf(n): return lpf_state(n, 2)

# ---------------------------------------------------------------- geometry (matches v2 builder)
NDIR=64
_ang=np.linspace(0,2*math.pi,NDIR,endpoint=False); _N=np.column_stack((np.cos(_ang),np.sin(_ang)))
def hp_arr(b): return np.sqrt(_N[:,0]**2+(b*_N[:,1])**2)              # container support
def hnew_arr(b,scale,ori):
    return scale*(np.sqrt(_N[:,0]**2+(b*_N[:,1])**2) if ori==0 else np.sqrt(_N[:,1]**2+(b*_N[:,0])**2))
def axes2(scale,ori,b):
    return (scale*scale,(scale*b)**2) if ori==0 else ((scale*b)**2,scale*scale)
_LAM=np.linspace(0,1,50)[1:-1]                                        # lambda grid (fine, matches certifier)
def pw_margin(dx,dy,ax_i,ax_j):
    cx=(1-_LAM)*ax_i[0]+_LAM*ax_j[0]; cy=(1-_LAM)*ax_i[1]+_LAM*ax_j[1]
    return np.max(_LAM*(1-_LAM)*(dx*dx/cx+dy*dy/cy))-1.0

# ---------------------------------------------------------------- spatial hash
class Grid:
    def __init__(self, cell): self.cell=cell; self.b={}
    def _k(self,x,y): return (int(math.floor(x/self.cell)), int(math.floor(y/self.cell)))
    def add(self,i,x,y): self.b.setdefault(self._k(x,y),[]).append(i)
    def near(self,x,y,reach):
        r=int(math.ceil(reach/self.cell)); k=self._k(x,y); out=[]
        for dx in range(-r,r+1):
            for dy in range(-r,r+1):
                out+=self.b.get((k[0]+dx,k[1]+dy),())
        return out

# ---------------------------------------------------------------- builder state
class State:
    BIG=0.05                                     # pieces with scale>=BIG are always tested
    def __init__(self, b, cn, ori, scales):
        self.b=b; self.scales=scales
        self.CN=list(map(list,cn)); self.ori=list(ori)         # normalized (x,u), orientations 0/1
        self.C=[(x, b*u) for (x,u) in cn]                      # actual centres
        self.ax=[axes2(scales[i],ori[i],b) for i in range(len(ori))]
        self.big=[i for i in range(len(ori)) if scales[i]>=self.BIG]
        self.grid=Grid(cell=max(2*self.BIG,0.04))
        for i in range(len(ori)):
            if scales[i]<self.BIG: self.grid.add(i,*self.C[i])
    def n(self): return len(self.ori)
    def margin(self, x, u, ori, k):
        """insertion margin for piece index k (scale scales[k]) at normalized (x,u), orientation ori."""
        b=self.b; cx,cy=x, b*u; s=self.scales[k]
        # containment (need centre inside; support test)
        hp=hp_arr(b); hn=hnew_arr(b,s,ori)
        proj=np.abs(cx*_N[:,0]+cy*_N[:,1])
        cont=float(np.min((hp-proj-hn)/hp))
        if cont<=0: return cont
        m=cont; axk=axes2(s,ori,b)
        cand=list(self.big)+self.grid.near(cx,cy,s+self.BIG)
        seen=set()
        for j in cand:
            if j in seen: continue
            seen.add(j)
            dx=cx-self.C[j][0]; dy=cy-self.C[j][1]
            if dx*dx+dy*dy > (s+self.scales[j])**2: continue     # cannot overlap
            pm=pw_margin(dx,dy,axk,self.ax[j])
            if pm<m: m=pm
            if m<=0: return m
        return m
    def true_margin(self, x, u, ori, k):
        """exact insertion margin against ALL placed pieces (vectorized) -- the value the
        certifier will see; used to gate acceptance so no pair is placed tighter than accept."""
        b=self.b; cx,cy=x,b*u; s=self.scales[k]
        hp=hp_arr(b); hn=hnew_arr(b,s,ori)
        cont=float(np.min((hp-np.abs(cx*_N[:,0]+cy*_N[:,1])-hn)/hp))
        if cont<=0 or self.n()==0: return cont
        axk=axes2(s,ori,b); Cn=np.asarray(self.C)
        dx=cx-Cn[:,0]; dy=cy-Cn[:,1]
        a0=np.array([a[0] for a in self.ax]); a1=np.array([a[1] for a in self.ax])
        L=_LAM[None,:]
        cxj=(1-L)*axk[0]+L*a0[:,None]; cyj=(1-L)*axk[1]+L*a1[:,None]
        pw=np.max(L*(1-L)*(dx[:,None]**2/cxj+dy[:,None]**2/cyj),axis=1)-1.0
        return min(cont, float(np.min(pw)))
    def truncate(self, m):
        return State(self.b, [tuple(c) for c in self.CN[:m]], self.ori[:m], self.scales)
    def prescreen(self, X, U, ori, k):
        """Vectorized lower-bounding margin over candidate arrays: containment plus PW against
        the (few) BIG pieces. Small-piece PW is added per-candidate afterwards, only for the
        top prescreen candidates, via .margin()."""
        b=self.b; s=self.scales[k]; hp=hp_arr(b); hn=hnew_arr(b,s,ori)
        CX=X; CY=b*U
        proj=np.abs(np.outer(CX,_N[:,0])+np.outer(CY,_N[:,1]))
        m=np.min((hp[None,:]-proj-hn[None,:])/hp[None,:],axis=1)          # containment
        if self.big:
            axk=axes2(s,ori,b); L=_LAM
            bx=np.array([self.C[j][0] for j in self.big]); by=np.array([self.C[j][1] for j in self.big])
            a0=np.array([self.ax[j][0] for j in self.big]); a1=np.array([self.ax[j][1] for j in self.big])
            dx=CX[:,None]-bx[None,:]; dy=CY[:,None]-by[None,:]           # Nc x nbig
            cxj=(1-L)[None,None,:]*axk[0]+L[None,None,:]*a0[None,:,None] # Nc x nbig x nlam
            cyj=(1-L)[None,None,:]*axk[1]+L[None,None,:]*a1[None,:,None]
            F=L[None,None,:]*(1-L)[None,None,:]*(dx[:,:,None]**2/cxj+dy[:,:,None]**2/cyj)
            pw=np.max(F,axis=2)-1.0                                       # Nc x nbig
            m=np.minimum(m, pw.min(axis=1))
        return m
    def place(self, x, u, ori, k):
        b=self.b; self.CN.append([x,u]); self.ori.append(ori)
        self.C.append((x,b*u)); self.ax.append(axes2(self.scales[k],ori,b))
        if self.scales[k]>=self.BIG: self.big.append(k)
        else: self.grid.add(k,x,b*u)
    def clone(self):
        st=object.__new__(State)
        st.b=self.b; st.scales=self.scales
        st.CN=[c[:] for c in self.CN]; st.ori=self.ori[:]
        st.C=self.C[:]; st.ax=self.ax[:]; st.big=self.big[:]
        g=Grid(self.grid.cell); g.b={key:v[:] for key,v in self.grid.b.items()}; st.grid=g
        return st

# ---------------------------------------------------------------- search one prime
def place_next(st, k, sob, accept, refine_iters, rng):
    # collect the best hash-margin candidates over both orientations
    pool=[]
    for ori in (0,1):
        b=st.b; s=st.scales[k]; hp=hp_arr(b); hn=hnew_arr(b,s,ori)
        X=sob[:,0]; U=sob[:,1]; CX=X; CY=b*U
        proj=np.abs(np.outer(CX,_N[:,0])+np.outer(CY,_N[:,1]))
        cont=np.min((hp[None,:]-proj-hn[None,:])/hp[None,:],axis=1)
        cand=np.where(cont>accept)[0]
        cand=cand[np.argsort(-cont[cand])][:1200]
        scored=[]
        for idx in cand:
            m=st.margin(X[idx],U[idx],ori,k)
            if m>accept: scored.append((m,X[idx],U[idx],ori))
        pool+=scored
    if not pool: return False
    pool.sort(key=lambda t:-t[0])
    # refine the top few on the (fast) hash margin, then GATE on the true margin
    best=(-9,None,None)
    for hm,x0,u0,ori in pool[:10]:
        f=lambda z:-st.margin(z[0],z[1],ori,k)
        r=minimize(f,[x0,u0],method='Nelder-Mead',
                   options={'maxiter':refine_iters,'xatol':1e-5,'fatol':1e-6})
        x,u=(r.x if -r.fun>hm else (x0,u0))
        tm=st.true_margin(float(x),float(u),ori,k)          # exact vs all placed
        if tm>best[0]: best=(tm,(float(x),float(u)),ori)
        if best[0]>max(accept*4,1e-3): break                # good enough, stop early
    if best[0]<=accept: return False
    (x,u),ori=best[1],best[2]; st.place(x,u,ori,k); return True

def top_placements(st, k, sob, accept, refine_iters, K=2, min_sep=0.12):
    """Up to K spatially-separated valid placements (true_margin, x, u, ori), best first.
    The beam keeps these non-greedy alternatives so a locally optimal but globally blocking
    placement does not doom the search."""
    pool=[]
    for ori in (0,1):
        X=sob[:,0]; U=sob[:,1]
        pm=st.prescreen(X,U,ori,k)                          # vectorized lower-bound margin
        good=np.where(pm>accept)[0]
        if good.size==0: continue
        good=good[np.argsort(-pm[good])][:400]              # refine only the top prescreen holes
        for idx in good:
            m=st.margin(X[idx],U[idx],ori,k)                # exact (adds small-piece PW)
            if m>accept: pool.append((m,X[idx],U[idx],ori))
    if not pool: return []
    pool.sort(key=lambda t:-t[0])
    picked=[]
    for hm,x0,u0,ori in pool:
        if len(picked)>=K: break
        if any((x0-px)**2+(u0-pu)**2<min_sep**2 for _,px,pu,po in picked if po==ori): continue
        f=lambda z:-st.margin(z[0],z[1],ori,k)
        r=minimize(f,[x0,u0],method='Nelder-Mead',
                   options={'maxiter':refine_iters,'xatol':1e-5,'fatol':1e-6})
        x,u=(r.x if -r.fun>hm else (x0,u0))
        tm=st.true_margin(float(x),float(u),ori,k)
        if tm>accept: picked.append((tm,float(x),float(u),ori))
    return picked

# ---------------------------------------------------------------- io
def sobol_field(count, seed):
    m=max(1,int(math.ceil(math.log2(max(2,count)))))
    return (2.0*qmc.Sobol(d=2,scramble=True,seed=seed).random_base2(m=m)-1.0)[:count]

def repair(st, k, n_floor, sob_count, seed_base, accept, refine_iters, rng, W=5, attempts=8):
    """Jam at piece k: back up the last W pieces (never below n_floor) and re-place the
    window k-W..k with fresh randomness, keeping the first success."""
    m=max(n_floor, k-W)
    for att in range(attempts):
        st2=st.truncate(m); ok=True
        for kk in range(m, k+1):
            sob=sobol_field(sob_count, seed_base + 9173*(att+1) + 31*kk)
            if not place_next(st2, kk, sob, accept, refine_iters, rng): ok=False; break
        if ok: return st2
    return None

def escalating_repair(st, k, n_floor, base_sob, seed, accept, refine_iters, rng, W0):
    """Try repair with geometrically larger backup windows until one escapes the jam."""
    for mult in (1, 2, 4, 8, 16, 32):
        W=W0*mult
        if k-W <= n_floor and mult>1 and W0*mult//2 >= (k-n_floor):
            pass  # already backing up to the floor; one more try then give up
        fixed=repair(st, k, n_floor, min(2*base_sob, 16384), seed+1000*mult,
                     accept, refine_iters, rng, W=W, attempts=8)
        if fixed is not None:
            return fixed, W
        if k-W <= n_floor:
            break
    return None, None

def load_start(path):
    d=json.load(open(path)); V=d.get('verification',d)
    if 'centers_normalized' in V:
        cn=np.asarray(V['centers_normalized'],float); b=float(V['aspect_ratio_b_over_a'])
        od=V.get('orientation_pattern_degrees'); ori=[1 if o==90 else 0 for o in od] if od else V['orientations']
    else:                                                    # checkpoint beam format
        st0=d['beam'][0]; cn=np.asarray(st0['centers_normalized'],float); b=float(st0['b']); ori=list(st0['orientations'])
    return b, cn, [int(o) for o in ori]

def save_head(path, st, primes, shares, state=2):
    b=st.b; C=[[x, b*u] for (x,u) in st.CN]; n=st.n()
    # area shares are exactly sigma_state(q) for the first n primes >= state; store as floats and
    # let the certifier reconstruct the exact rationals from (state, n) -- their exact denominators
    # run to thousands of digits and must never be str()'d.
    obj={"format":"ellipse-binary-head","status":"BUILT (verify with ellipse_certify.py)",
         "state":int(state),"n":n,"q_last":int(primes[n-1]),
         "verification":{
             "state":int(state),
             "aspect_ratio_b_over_a":b,
             "orientation_pattern_degrees":[90 if o else 0 for o in st.ori],
             "centers":C,"centers_normalized":st.CN,
             "area_shares_float":[float(shares[i]) for i in range(n)],
             "sum_head_area_float":sum(float(shares[i]) for i in range(n))}}
    tmp=path+".tmp"
    with open(tmp,"w") as f: json.dump(obj,f)
    os.replace(tmp,path)                                    # atomic: never leaves a partial file

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resume",default=None); ap.add_argument("--target-n",type=int,default=1650)
    ap.add_argument("--state",type=int,default=2,help="least admissible child prime; 2 = container")
    ap.add_argument("--rho",type=float,default=0.6296271949192772,help="aspect ratio (fresh start)")
    ap.add_argument("--out",default="ellipse_head.json"); ap.add_argument("--checkpoint",default="ellipse_ckpt.json")
    ap.add_argument("--sobol",type=int,default=3072); ap.add_argument("--seed",type=int,default=7)
    ap.add_argument("--accept",type=float,default=2e-5); ap.add_argument("--refine-iters",type=int,default=120)
    ap.add_argument("--checkpoint-every",type=int,default=25)
    ap.add_argument("--repair-window",type=int,default=5)
    ap.add_argument("--beam",type=int,default=6)
    ap.add_argument("--per-parent",type=int,default=2)
    a=ap.parse_args()

    primes,shares,scales=lpf_state(a.target_n, a.state)
    if len(scales)<a.target_n: raise SystemExit("need more primes")
    if a.resume:
        b,cn,ori=load_start(a.resume); start=State(b,cn,ori,scales)
    else:                                                    # fresh: beam places everything, incl.
        b=a.rho; start=State(b, [], [], scales)              # the spread-out big pieces, from empty
    n_floor=start.n()
    st_p = f" (through p={primes[start.n()-1]})" if start.n()>0 else ""
    print(f"state {a.state}: start n={start.n()}{st_p}; target {a.target_n}; "
          f"b={b:.6f}; beam={a.beam}, per-parent={a.per_parent}")
    t0=time.time()
    beam=[start]; stalls=0
    while beam[0].n()<a.target_n:
        k=beam[0].n()
        big = scales[k] >= 0.12                                # big piece: search wide for offsets
        K = max(a.per_parent, 8) if big else a.per_parent
        msep = 0.40 if big else 0.12
        sob=sobol_field(a.sobol, a.seed+k+stalls*100003)      # reshuffle field on a stall
        children=[]                                           # (tm, parent, x, u, ori)
        for state in beam:
            for (tm,x,u,ori) in top_placements(state,k,sob,a.accept,a.refine_iters,
                                               K=K,min_sep=msep):
                children.append((tm,state,x,u,ori))
        if not children:
            stalls+=1
            if stalls<=10: continue                           # try fresh fields before giving up
            print(f"  JAM at piece {k} (prime {primes[k]}): beam exhausted after reshuffles. "
                  f"Resume from checkpoint with a different --seed, or raise --beam/--sobol; "
                  f"checkpoint saved.")
            save_head(a.checkpoint,beam[0],primes,shares,a.state); break
        stalls=0
        children.sort(key=lambda t:-t[0])
        beam=[p.clone() for (_,p,_,_,_) in children[:a.beam]]
        for st,(tm,parent,x,u,ori) in zip(beam, children[:a.beam]):
            st.place(x,u,ori,k)
        if beam[0].n()%a.checkpoint_every==0 or beam[0].n()==a.target_n:
            save_head(a.checkpoint,beam[0],primes,shares,a.state)
            print(f"  n={beam[0].n():>4} (p={primes[beam[0].n()-1]:>6})  {time.time()-t0:6.1f}s  "
                  f"beam={len(beam)}  -> checkpoint")
    save_head(a.out,beam[0],primes,shares,a.state)
    print(f"wrote {a.out} at n={beam[0].n()} ({time.time()-t0:.1f}s). "
          f"Verify: python3 ellipse_certify.py certify --in {a.out}")

if __name__=="__main__":
    main()
