#!/usr/bin/env python3
"""Exploratory axis-parallel rectangle aspect-ratio sweep (NOT a certificate).

Run without arguments for an overview, or use --count 5000 --steps 101.
This tests explicit finite guillotine packing heuristics, and reports failures;
failure of a heuristic is not an impossibility result.
"""
import argparse
import csv
import math


def primes_through_count(count):
    cap = max(100, int(2.2*count*math.log(max(count, 4))))
    while True:
        sieve = bytearray(b'\x01') * (cap+1)
        sieve[:2] = b'\x00\x00'
        for p in range(2, math.isqrt(cap)+1):
            if sieve[p]:
                sieve[p*p:cap+1:p] = b'\x00' * ((cap-p*p)//p+1)
        ps = [p for p in range(2,cap+1) if sieve[p]]
        if len(ps)>=count: return ps[:count]
        cap*=2


def attempt(r, primes, method):
    W, H = math.sqrt(r), 1/math.sqrt(r)
    free = [(W,H)]
    R = 1.0
    for j,p in enumerate(primes):
        d = R/p
        a,b = math.sqrt(r*d), math.sqrt(d/r)
        choices = []
        for i,(w,h) in enumerate(free):
            for k,(x,y) in enumerate(((a,b),(b,a))):
                if w+1e-13<x or h+1e-13<y:continue
                slackx,slacky=max(0,w-x),max(0,h-y)
                if method=='tight':
                    key=(min(slackx,slacky),w*h, max(slackx,slacky))
                elif method=='area':
                    key=(w*h, min(slackx,slacky), max(slackx,slacky))
                else:
                    key=(max(slackx,slacky),w*h,min(slackx,slacky))
                choices.append((key,i,k))
        if not choices: return (j,p,len(free),R)
        _,i,k=min(choices)
        w,h=free.pop(i)
        x,y=((a,b),(b,a))[k]
        # Disjoint guillotine residuals: right strip and upper strip.
        if w-x>1e-14 and h>1e-14:free.append((w-x,h))
        if h-y>1e-14 and x>1e-14:free.append((x,h-y))
        R*=(p-1)/p
    return (len(primes),None,len(free),R)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--count',type=int,default=1200)
    ap.add_argument('--steps',type=int,default=81)
    ap.add_argument('--min-r',type=float,default=1/math.sqrt(2)+1/math.sqrt(6))
    ap.add_argument('--max-r',type=float,default=math.sqrt(6))
    ap.add_argument('--csv',default='rectangle_aspect_sweep.csv')
    args=ap.parse_args()
    ps=primes_through_count(args.count)
    rows=[]
    for k in range(args.steps):
        r=args.min_r+(args.max_r-args.min_r)*k/(args.steps-1)
        outcomes=[(m,*attempt(r,ps,m)) for m in ('tight','area','maxslack')]
        best=max(outcomes,key=lambda z:z[1])
        rows.append({'r':f'{r:.10f}','packed':best[1],
                     'next_prime':best[2] or '', 'strategy':best[0],
                     'free_rectangles':best[3], 'remaining_area':best[4]})
        print(f'{k+1:3d}/{args.steps}: r={r:.7f}  {best[1]}/{args.count}  '
              f'best={best[0]}  next={best[2]}',flush=True)
    with open(args.csv,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('saved',args.csv)


if __name__=='__main__':main()
