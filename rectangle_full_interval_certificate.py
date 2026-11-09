#!/usr/bin/env python3
"""Finite-cover certificate for the COMPLETE allowable aspect-ratio interval.

Run without flags. Requires mpmath and the local rectangle_aspect helpers.
All 16 head schedules pack the first 1282 primes (through p=10477);
every piece after p=3 is aligned with the container. For each schedule,
160-bit interval inequalities certify an entire r-interval, not one r.
Exact endpoint contact identities are handled algebraically.
Every post-head prime through 10^7 is evaluated in outward-rounded interval
arithmetic. The infinite range uses the Appendix-A rho bound and explicit
Rosser--Schoenfeld/Dusart estimates at 10^7.
"""
import json
import math
import subprocess
import sys
from decimal import Decimal
from mpmath import iv
from rectangle_aspect_path_range import path_range,certify_exact_endpoint
from rectangle_aspect_sweep import primes_through_count


COUNT=1282
LAST_HEAD_PRIME=10477
X1=10_000_000
SCHEDULES=[
    ('1.11841332','RR'),
    ('1.118418','RR'),
    ('1.11856','RR'),
    ('1.11887','RR'),
    ('1.1189','RR'),
    ('1.119','RR'),
    ('1.18','RR'),
    ('1.2','RR'),
    ('1.3','RA'),
    ('1.38','RA'),
    ('1.4','RA'),
    ('1.65','AR'),
    ('2.0','AR'),
    ('2.02','AR'),
    ('2.4','AR'),
    ('2.44948','AR'),
]


def primes_up_to(n):
    """Exact integer sieve; no prime-density approximation in the finite window."""
    sieve=bytearray(b'\x01')*(n+1)
    sieve[0:2]=b'\x00\x00'
    for p in range(2,math.isqrt(n)+1):
        if sieve[p]:
            sieve[p*p:n+1:p]=b'\x00'*((n-p*p)//p+1)
    return (p for p in range(2,n+1) if sieve[p])


def certify_uniform_mixed_tail():
    iv.prec=160
    half=1/iv.sqrt(iv.mpf(2))
    third=1/iv.sqrt(iv.mpf(6))
    cmax=iv.sqrt(iv.mpf(6))+1/iv.sqrt(iv.mpf(6))
    excess=(cmax-2)*(half+third)
    R,T=iv.mpf(1),iv.mpf(0)
    maximum=(-1,None)
    maximum_late=(-1,None)
    checked=0
    prime_total=0
    for j,p in enumerate(primes_up_to(X1)):
        if j==COUNT-1 and p!=LAST_HEAD_PRIME:
            raise AssertionError('finite-head prime cutoff does not match')
        if p>LAST_HEAD_PRIME:
            B=(2+2*T+excess)/iv.sqrt(p*R)+iv.mpf(j-1)/p
            if not (B.b<iv.mpf(1).a):
                raise AssertionError(f'uniform mixed tail failed at {p}: {B}')
            if float(B.b)>maximum[0]:maximum=(float(B.b),p)
            if p>50021 and float(B.b)>maximum_late[0]:
                maximum_late=(float(B.b),p)
            checked+=1
        T+=iv.sqrt(R/p)
        R*=iv.mpf(p-1)/p
        prime_total=j+1
    if prime_total!=664579:
        raise AssertionError(f'unexpected pi(10^7)={prime_total}')
    # For q>=10^7, Appendix A proves rho(q)=(3*T_inclusive+sqrt(2))
    # /sqrt(q*R_<q)<0.614. RS gives R_<q>=M(q-1)>=beta/log(q)
    # with conservative beta evaluated at Q1=10^5. The latter two
    # decreasing bounds are evaluated at X1 in outward-rounded arithmetic.
    L=iv.log(iv.mpf(X1))
    L0=iv.log(iv.mpf(100000))
    beta=iv.exp(-iv.euler)*(1-1/L0**2)
    inv_sqrt=iv.sqrt(L/(beta*X1))
    count_bound=(1+iv.mpf('1.2762')/L)/L
    analytic=(iv.mpf(2)/3)*iv.mpf('0.614')+(
        2+excess-(iv.mpf(2)/3)*iv.sqrt(iv.mpf(2))
    )*inv_sqrt+count_bound
    if not (analytic.b<iv.mpf(1).a):
        raise AssertionError(f'infinite analytic majorant failed: {analytic}')
    return {'finite_primes_checked':checked,'max_finite_upper':maximum[0],
            'max_finite_at_prime':maximum[1],
            'post_50021_max_upper':maximum_late[0],
            'post_50021_max_at_prime':maximum_late[1],
            'window_endpoint':X1,'all_primes_in_window':prime_total,
            'beta_lower':float(beta.a),
            'inverse_sqrt_upper_at_X1':float(inv_sqrt.b),
            'prime_count_upper_at_X1':float(count_bound.b),
            'infinite_majorant_upper':float(analytic.b)}


def main():
    iv.prec=200
    intervals=[]
    for index,(r,pair) in enumerate(SCHEDULES):
        lowroot,highroot,(lo,hi),signature=path_range(
            float(r),COUNT,prec=200,pair=pair,method='area',return_signature=True)
        if index==0:
            margin,n,endpoint_sig=certify_exact_endpoint(float(r),pair,'lower',COUNT)
            if endpoint_sig!=signature:raise AssertionError('lower endpoint schedule mismatch')
            print(f'LOWER ENDPOINT exact-contact check: {n} strict inequalities, min={margin:.3g}')
        if index==len(SCHEDULES)-1:
            margin,n,endpoint_sig=certify_exact_endpoint(float(r),pair,'upper',COUNT)
            if endpoint_sig!=signature:raise AssertionError('upper endpoint schedule mismatch')
            print(f'UPPER ENDPOINT exact-contact check: {n} strict inequalities, min={margin:.3g}')
        intervals.append({'sample':r,'pair':pair,'lo':repr(lo),'hi':repr(hi),
                          'sha256_choices':signature,'pieces':COUNT})
        print(f'PATH {index+1:02d}/{len(SCHEDULES)} r={r} {pair}: [{lo:.15f}, {hi:.15f}]',flush=True)
    # Decimal endpoints are exact rational bounds already interval-certified.
    sorted_intervals=sorted(intervals,key=lambda item:Decimal(item['lo']))
    reached=Decimal(sorted_intervals[0]['hi'])
    for item in sorted_intervals[1:]:
        left,right=Decimal(item['lo']),Decimal(item['hi'])
        if left>reached:
            raise AssertionError(f'UNSUPPORTED ASPECT GAP ({reached},{left})')
        reached=max(reached,right)
    # The first/last path extend to the algebraic endpoints via their single
    # symbolically nonnegative contact inequalities; all other constraints
    # were interval-certified strictly positive at the exact endpoints.
    first=intervals[0]
    last=intervals[-1]
    if sorted_intervals[0]['sample']!=first['sample']:
        raise AssertionError('lower-endpoint schedule not first')
    if reached!=Decimal(last['hi']):
        raise AssertionError('upper-endpoint schedule not last')
    print('CLOSED-INTERVAL COVER: CERTIFIED (16 schedules plus two exact contacts)')

    tail=certify_uniform_mixed_tail()
    print('UNIFORM MIXED-ORIENTATION TAIL: CERTIFIED',tail)
    proc=subprocess.run([sys.executable,'rectangle_aspect_recursive_certificate.py',
                         '--r','1','--square-states'],capture_output=True,text=True)
    if proc.returncode or 'OVERALL: CERTIFIED finite recursive inputs' not in proc.stdout:
        raise AssertionError('universal square-state family failed:\n'+proc.stdout+proc.stderr)
    print('ALL SHIFTED STATES s>=3: CERTIFIED (universal aligned-square family)')
    record={'status':'CERTIFIED complete interval','head_pieces':COUNT,
            'head_last_prime':LAST_HEAD_PRIME,'schedules':intervals,
            'tail':tail,'universal_states':'CERTIFIED',
            'aspect_interval':['1/sqrt(2)+1/sqrt(6)','sqrt(6)']}
    with open('rectangle_full_interval_certificate.json','w') as f:
        json.dump(record,f,indent=2)
    print('OVERALL: CERTIFIED COMPLETE ASPECT INTERVAL')
    print('Saved rectangle_full_interval_certificate.json')


if __name__=='__main__':main()
