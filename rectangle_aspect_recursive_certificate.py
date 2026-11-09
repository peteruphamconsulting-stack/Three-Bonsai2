#!/usr/bin/env python3
"""Verify all required shifted-state finite heads at one fixed aspect ratio.

Run: python3 rectangle_aspect_recursive_certificate.py --r 1.2
This is a preliminary deterministic interval certificate. It uses the common
analytic continuation bound after q=50000, valid on the two-piece r interval.
The float crossover search is followed by fail-closed interval checks.
"""
import argparse
import math
from mpmath import iv
from rectangle_aspect_head_certificate import verify,aspect_value
from rectangle_aspect_sweep import primes_through_count


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--r',required=True)
    ap.add_argument('--prec',type=int,default=160)
    ap.add_argument('--square-states',action='store_true',
                    help='certify universal aligned-only shifted states s>=3 at r=1')
    args=ap.parse_args()
    iv.prec=args.prec
    rir,rv=aspect_value(args.r)
    if not args.square_states and not (1/math.sqrt(2)+1/math.sqrt(6)-1e-12<=rv<=math.sqrt(6)+1e-12):
        raise SystemExit('Aspect is outside the two-piece feasible interval.')
    P=primes_through_count(5200)
    P=[p for p in P if p<=50021]
    R=[1.0];T=[0.0]
    for p in P:
        d=R[-1]/p
        T.append(T[-1]+math.sqrt(d))
        R.append(R[-1]*(p-1)/p)
    C=rv+1/rv
    state_list=[]
    for i,s in enumerate(P):
        if args.square_states and i==0:continue
        bmax=-1.; lastbad=-1
        for j in range(i,len(P)):
            p=P[j]
            B=(2*math.sqrt(R[i])+C*(T[j]-T[i]))/math.sqrt(p*R[j])+(j-i-1)/p
            if B>=1-1e-9:lastbad=j
            bmax=max(B,bmax)
        if lastbad>=0:state_list.append((i,lastbad,bmax))
        elif i>100:break
    print(f'aspect {rv:.10f}: {len(state_list)} state heads needed (float search), '
          f'last state={P[state_list[-1][0]] if state_list else None}',flush=True)
    errors=0
    for i,last,bmax in state_list:
        s=P[i]
        count=last-i+1
        result=verify(args.r,count,args.prec,state=s)
        if not result[0]:
            print(f'FAILED HEAD state {s} at piece {result[2]}',flush=True)
            errors+=1
            continue
        # Recompute the normalized criterion directly in interval arithmetic.
        R_iv=iv.mpf(1)
        T_iv=iv.mpf(0)
        Rs=Ts=None
        total_count=0
        uppermax=-1.
        for j,p in enumerate(P):
            if j==i:Rs,Ts=R_iv,T_iv
            if j>last and Rs is not None:
                br=(2*iv.sqrt(Rs)+(rir+1/rir)*(T_iv-Ts))/iv.sqrt(p*R_iv)+iv.mpf(j-i-1)/p
                if not (br.b<iv.mpf(1).a):
                    errors+=1
                    print(f'FAILED CAPACITY state {s} at p={p}, upper={float(br.b)}',flush=True)
                    break
                uppermax=max(uppermax,float(br.b));total_count+=1
            T_iv+=iv.sqrt(R_iv/p)
            R_iv*=iv.mpf(p-1)/p
        else:
            print(f'CERTIFIED state {s}: {count} pieces through {P[last]}, '
                  f'min head margin={result[3]:.3g}, tail count={total_count}, '
                  f'capacity max={uppermax:.9g}',flush=True)
    # Certify the first state needing no head. The monotone comparison then
    # transfers this check to every later state, without trusting a float scan.
    terminal_index=state_list[-1][0]+1 if state_list else 0
    terminal=P[terminal_index]
    R_iv=iv.mpf(1);T_iv=iv.mpf(0)
    Rs=Ts=None
    terminal_max=-1.; terminal_checks=0
    for j,p in enumerate(P):
        if j==terminal_index:Rs,Ts=R_iv,T_iv
        if j>=terminal_index:
            B=(2*iv.sqrt(Rs)+(rir+1/rir)*(T_iv-Ts))/iv.sqrt(p*R_iv)+iv.mpf(j-terminal_index-1)/p
            if not (B.b<iv.mpf(1).a):
                errors+=1
                print(f'FAILED terminal head-free state {terminal}, prime {p}',flush=True)
                break
            terminal_max=max(terminal_max,float(B.b));terminal_checks+=1
        T_iv+=iv.sqrt(R_iv/p)
        R_iv*=iv.mpf(p-1)/p
    else:
        print(f'CERTIFIED all states >= {terminal}: head-free state checked at '
              f'{terminal_checks} primes, max upper={terminal_max:.9g}',flush=True)
    print('OVERALL:', 'CERTIFIED finite recursive inputs' if errors==0 else f'{errors} failures')
    raise SystemExit(0 if errors==0 else 1)


if __name__=='__main__':main()
