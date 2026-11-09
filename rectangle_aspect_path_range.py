#!/usr/bin/env python3
"""Determine the exact-in-r feasibility interval of one fixed guillotine path.

The path is chosen at --r with floats; coefficients of all fit inequalities
are reconstructed in 160-bit interval arithmetic and the reported inner
interval is checked at both endpoints (linear inequalities in r).
"""
import argparse
import hashlib
import math
from mpmath import iv
from rectangle_aspect_sweep import primes_through_count


def mid(x):return float(x.mid)


def path_range(rv,count=5133,prec=160,pair=None,method='tight',return_signature=False):
    iv.prec=prec
    # normalized dimensions: x=A+B/r, y=C+D*r.
    free=[(iv.mpf(1),iv.mpf(0),iv.mpf(1),iv.mpf(0))]
    shadows=[(1.,1.)]
    R=iv.mpf(1)
    constraints=[]
    choices=[]
    for j,p in enumerate(primes_through_count(count)):
        d=R/p; t=iv.sqrt(d); tf=mid(t)
        possible=[]
        for i,(wf,hf) in enumerate(shadows):
            allowed=((1 if pair[0]=='R' else 0),) if pair and p==2 else (((1 if pair[1]=='R' else 0),) if pair and p==3 else ((0,) if pair else (0,1)))
            for k in allowed:
                xf,yf=(tf,tf) if k==0 else (tf/rv,tf*rv)
                if wf>=xf-1e-12 and hf>=yf-1e-12:
                    if method=='tight':key=(min(wf-xf,hf-yf),wf*hf,max(wf-xf,hf-yf))
                    elif method=='area':key=(wf*hf,min(wf-xf,hf-yf),max(wf-xf,hf-yf))
                    else:key=(max(wf-xf,hf-yf),wf*hf,min(wf-xf,hf-yf))
                    possible.append((key,i,k))
        if not possible:raise RuntimeError(f'failed at p={p}, j={j}')
        _,i,k=min(possible)
        choices.append((p,i,k))
        A,B,C,D=free.pop(i)
        wf,hf=shadows.pop(i)
        a,b,c,e=(t,iv.mpf(0),t,iv.mpf(0)) if k==0 else (iv.mpf(0),t,iv.mpf(0),t)
        xf,yf=(tf,tf) if k==0 else (tf/rv,tf*rv)
        constraints.extend(((p,'x',A-a,B-b),(p,'y',D-e,C-c)))
        # right strip: (w-x,h), upper-left: (x,h-y)
        free.append((A-a,B-b,C,D))
        free.append((a,b,C-c,D-e))
        shadows.append((wf-xf,hf))
        shadows.append((xf,hf-yf))
        R*=iv.mpf(p-1)/p
    low,high=0.,float('inf')
    for _,_,a,b in constraints:
        af,bf=mid(a),mid(b)
        if abs(af)<1e-25:
            if bf< -1e-10:raise RuntimeError('inconsistent constant inequality')
            continue
        root=-bf/af
        if af>0:low=max(low,root)
        else:high=min(high,root)
    # Find a conservative interval inside the true feasible segment.
    eps=1e-9*max(1.,abs(low),abs(high) if math.isfinite(high) else 1.)
    lo=max(1.,low+eps)
    hi=min(3.,high-eps)
    if lo>hi:return low,high,None
    for endpoint in (lo,hi):
        v=iv.mpf(repr(endpoint))
        for _,_,a,b in constraints:
            if not (a*v+b).a>=0:
                raise RuntimeError(f'endpoint interval failed at {endpoint} for {mid(a)}r+{mid(b)}')
    result=(low,high,(lo,hi))
    if return_signature:
        return (*result,hashlib.sha256(repr(choices).encode()).hexdigest())
    return result


def certify_exact_endpoint(rv,pair,which,count=1282,prec=200):
    """Check an endpoint-specific schedule; skip only its known exact p3 contact."""
    iv.prec=prec
    exact=(1/iv.sqrt(iv.mpf(2))+1/iv.sqrt(iv.mpf(6))) if which=='lower' else iv.sqrt(iv.mpf(6))
    free=[(iv.mpf(1),iv.mpf(0),iv.mpf(1),iv.mpf(0))]
    R=iv.mpf(1)
    constraints=[]
    choices=[]
    shadows=[(1.,1.)]
    for p in primes_through_count(count):
        t=iv.sqrt(R/p);tf=mid(t)
        allowed=((1 if pair[0]=='R' else 0),) if p==2 else (((1 if pair[1]=='R' else 0),) if p==3 else (0,))
        candidates=[]
        for i,(wf,hf) in enumerate(shadows):
            for k in allowed:
                xf,yf=(tf,tf) if k==0 else (tf/rv,tf*rv)
                if wf>=xf-1e-12 and hf>=yf-1e-12:
                    candidates.append(((wf*hf,min(wf-xf,hf-yf),max(wf-xf,hf-yf)),i,k))
        if not candidates:raise RuntimeError(f'failed candidate at p={p}')
        _,i,k=min(candidates)
        choices.append((p,i,k))
        A,B,C,D=free.pop(i);wf,hf=shadows.pop(i)
        a,b,c,e=(t,iv.mpf(0),t,iv.mpf(0)) if k==0 else (iv.mpf(0),t,iv.mpf(0),t)
        xf,yf=(tf,tf) if k==0 else (tf/rv,tf*rv)
        constraints.extend(((p,'x',A-a,B-b),(p,'y',D-e,C-c)))
        free.extend(((A-a,B-b,C,D),(a,b,C-c,D-e)))
        shadows.extend(((wf-xf,hf),(xf,hf-yf)))
        R*=iv.mpf(p-1)/p
    contact_axis='x' if which=='lower' else 'y'
    expected_first=((2,0,1),(3,0,1)) if which=='lower' else ((2,0,0),(3,0,1))
    if tuple(choices[:2])!=expected_first:
        raise RuntimeError(f'endpoint contact not in expected right-strip configuration: {choices[:2]}')
    # In the lower RR configuration, the p3 x-slack is
    #   1 - (1/sqrt(2)+1/sqrt(6))/r, exactly zero at r=r_min.
    # In the upper AR configuration, the p3 y-slack is
    #   1-r/sqrt(6), exactly zero at r=sqrt(6).
    strict_min=float('inf')
    for p,axis,a,b in constraints:
        if p==3 and axis==contact_axis:continue
        val=a*exact+b
        if not val.a>0:
            raise RuntimeError(f'non-contact inequality at p={p} axis={axis}: {val}')
        strict_min=min(strict_min,float(val.a))
    return strict_min,len(constraints)-1,hashlib.sha256(repr(choices).encode()).hexdigest()


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r',type=float,required=True)
    ap.add_argument('--count',type=int,default=5133)
    ap.add_argument('--pair',choices=['RR','RA','AR'],default=None,
                    help='force first two orientations and align all later pieces')
    ap.add_argument('--method',choices=['tight','area','maxslack'],default='tight')
    a=ap.parse_args()
    print('sample',a.r,'pair',a.pair,'method',a.method,'path interval',path_range(a.r,a.count,pair=a.pair,method=a.method),flush=True)
