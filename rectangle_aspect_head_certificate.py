#!/usr/bin/env python3
"""Interval-verify a finite guillotine head for an arbitrary fixed aspect r.

This certifies only the state-2 head; the recursive shifted-state family needs
separate heads unless r=sqrt(2). Use --r 1.2 (decimal); default count=5133.
Heuristic choices use float; all actual geometric acceptance uses mpmath.iv.
"""
import argparse
import math
from mpmath import iv
from rectangle_aspect_sweep import primes_through_count


def upper(x): return float(x.b)
def lower(x): return float(x.a)


def aspect_value(value):
    if value=='sqrt6':return iv.sqrt(iv.mpf(6)),math.sqrt(6)
    if value=='lower':
        return 1/iv.sqrt(iv.mpf(2))+1/iv.sqrt(iv.mpf(6)),1/math.sqrt(2)+1/math.sqrt(6)
    r=iv.mpf(value)
    return r,float(r.mid)


def verify(r_string, count, prec=160, state=2):
    iv.prec=prec
    r,rf=aspect_value(r_string)
    W,H=iv.sqrt(r),1/iv.sqrt(r)
    fw,fh=math.sqrt(rf),1/math.sqrt(rf)
    free=[(W,H)]
    shadow=[(fw,fh)]
    rem=iv.mpf(1)
    min_gap=float('inf')
    if state==2 and r_string in ('sqrt6','lower'):
        # The two-piece necessary-condition endpoints have an exact contact.
        # Specify their disjoint residual strips algebraically, rather than
        # attempting to certify a zero interval difference as positive.
        half=1/iv.sqrt(iv.mpf(2))
        third=1/iv.sqrt(iv.mpf(6))
        if r_string=='lower':
            # Both pieces rotated, side by side across the whole width.
            free=[(H*half,H-W*half),(H*third,H-W*third)]
        else:
            # Piece 2 aligned at left-bottom; piece 3 rotated at right,
            # spanning the full container height exactly.
            free=[(W*half,H*(1-half)),(W*(1-half)-H*third,H)]
        if any(x.a<=0 or y.a<=0 for x,y in free):
            raise RuntimeError('symbolic endpoint residual is not positive')
        shadow=[(float(x.mid),float(y.mid)) for x,y in free]
        rem=iv.mpf(1)/3
    all_primes=primes_through_count(count+max(5,2*state))
    state_primes=[p for p in all_primes if p>=state][:count]
    if len(state_primes)!=count:raise ValueError('prime sieve did not include full state head')
    for n,p in enumerate(state_primes,1):
        if state==2 and r_string in ('sqrt6','lower') and p in (2,3):
            continue
        d=rem/p
        a,b=iv.sqrt(r*d),iv.sqrt(d/r)
        af,bf=math.sqrt(rf*float(d.mid)),math.sqrt(float(d.mid)/rf)
        candidates=[]
        for i,(w,h) in enumerate(shadow):
            for k,(x,y) in enumerate(((af,bf),(bf,af))):
                if w>=x-1e-12 and h>=y-1e-12:
                    candidates.append(((min(w-x,h-y),w*h,max(w-x,h-y)),i,k))
        success=False
        for _,i,k in sorted(candidates):
            w,h=free[i]
            x,y=(a,b) if k==0 else (b,a)
            if not (w.a>=x.b and h.a>=y.b):continue
            min_gap=min(min_gap,lower(w-x),lower(h-y))
            wf,hf=shadow[i]
            xf,yf=(af,bf) if k==0 else (bf,af)
            free.pop(i);shadow.pop(i)
            for dims,fl in [((w-x,h),(wf-xf,hf)),((x,h-y),(xf,hf-yf))]:
                if dims[0].a>0 and dims[1].a>0:
                    free.append(dims);shadow.append(fl)
            success=True
            break
        if not success:
            return False,n,p,min_gap,len(free)
        rem*=iv.mpf(p-1)/p
    ledger=sum((w*h for w,h in free),iv.mpf(0))
    if not (ledger.a<=rem.b and rem.a<=ledger.b):
        raise RuntimeError('Free-area ledger not enclosed')
    return True,count,None,min_gap,len(free)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r',required=True,help='decimal ratio r >= 1')
    ap.add_argument('--count',type=int,default=5133)
    ap.add_argument('--prec',type=int,default=160)
    args=ap.parse_args()
    ok,n,p,gap,nfree=verify(args.r,args.count,args.prec)
    print(f'r={args.r}, status={"CERTIFIED FINITE HEAD" if ok else "FAILED HEAD"}, '
          f'pieces={n if ok else n-1}/{args.count}, first_failure={p}, '
          f'min_margin_lower={gap:.9g}, free_rectangles={nfree}')
    raise SystemExit(0 if ok else 1)


if __name__=='__main__':main()
