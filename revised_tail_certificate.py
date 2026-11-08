#!/usr/bin/env python3
"""Outward-rounded LPF packing tail certificates; no geometric head regeneration.

Usage:
  python3 revised_tail_certificate.py ellipse
  python3 revised_tail_certificate.py triangle
  python3 revised_tail_certificate.py triangle --design uniform8 --X1 10000000

The finite prime window is checked in mpmath.iv arithmetic.  The infinite tail
uses the Rosser--Schoenfeld/Dusart majorant with its *full* prefix term and
prime-count term.  FAIL is a real failure of this sufficient bound, not a
counterexample to the existence of a packing.
"""
import argparse
import math
from mpmath import iv


def primes_up_to(n):
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[0:2] = b'\x00\x00'
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p:n + 1:p] = b'\x00' * ((n - p * p) // p + 1)
    return (p for p in range(2, n + 1) if sieve[p])


def low(x):
    return float(x.a)


def high(x):
    return float(x.b)


def ellipse_coeff(rho):
    """Rigorous enclosing series for E(m), m=1-rho**4.

    E(m)=(pi/2)*sum_k [(-1)^k binom(1/2,k) binom(2k,k)/4^k] m^k.
    All terms after k=0 are negative; their magnitudes are bounded by m^k.
    Thus adding [-m^(K+1)/(1-m),0] to the partial sum encloses the tail.
    """
    r = iv.mpf(str(rho))
    m = 1 - r**4
    if not (0 < low(m) <= high(m) < 1):
        raise ValueError('ellipse aspect must lie strictly between zero and one')
    summand = iv.mpf(1)
    total = summand
    K = 240
    for k in range(1, K + 1):
        # coefficient recurrence: c_k/c_{k-1} = ((k-1/2)(k-1/2))/(k*k)
        # with c_1=-1/4 and c_k negative for all k>=1.
        if k == 1:
            summand = -m / 4
        else:
            summand *= m * iv.mpf(2*k - 3) * iv.mpf(2*k - 1) / (4*k*k)
        total += summand
    tail_iv = m**(K+1) / (1-m)
    E_low = (iv.pi/2) * (total - tail_iv)
    E_high = (iv.pi/2) * total
    E = iv.mpf([E_low.a, E_high.b])
    return 2*E/(iv.pi*r)


def certify(model='ellipse', design='disciplined', rho='0.6296271949192772',
            p0=1000, Q1=100000, X1=10000000, prec=128,
            head_prime=None, quiet=False):
    if X1 <= Q1 or Q1 < 100000:
        raise ValueError('Require X1>Q1>=100000 for the cited explicit bounds')
    if model == 'ellipse' and design != 'disciplined':
        raise ValueError('ellipse uses its two-orientation averaged criterion')
    iv.prec = prec
    R = iv.mpf(1)
    Tstrict = iv.mpf(0)
    Cf = iv.mpf(0)
    N = 0
    if model == 'ellipse':
        kstar = ellipse_coeff(rho)
        coeff = 1 + kstar  # 2*kbar
        head_prime = 13967 if head_prime is None else head_prime
        print('ellipse kappa*: [%.12f, %.12f]; kbar: [%.12f, %.12f]' %
              (low(kstar), high(kstar), low(coeff/2), high(coeff/2)))
    elif model == 'triangle':
        head_prime = (115363 if design == 'disciplined' else 662177) if head_prime is None else head_prime
        coeff = iv.mpf(3) if design == 'disciplined' else 3*iv.sqrt(2)
    else:
        raise ValueError('Unknown model')

    max_after = (-1.0, None)
    last_fail = None
    first_after_head = None
    T_Q1 = None
    Cf_Q1 = None
    for p in primes_up_to(X1):
        root_d = iv.sqrt(R / p)
        if model == 'ellipse':
            numerator = coeff * (1 + Tstrict)
        elif design == 'disciplined':
            numerator = 3 + 3*Tstrict + Cf
        else:
            numerator = 2 + coeff*Tstrict
        criterion = numerator / iv.sqrt(p*R) + iv.mpf(N)/p
        if p > head_prime:
            if first_after_head is None:
                first_after_head = (p, high(criterion))
            upper = high(criterion)
            if upper > max_after[0]:
                max_after = (upper, p)
            if not (criterion.b < iv.mpf(1).a):
                last_fail = (p, low(criterion), high(criterion))
        if p <= Q1:
            T_Q1 = Tstrict + root_d
            if model == 'triangle' and design == 'disciplined' and p <= p0:
                Cf_Q1 = Cf + root_d
        Tstrict += root_d
        if model == 'triangle' and design == 'disciplined' and p <= p0:
            Cf += root_d
        R *= iv.mpf(p - 1) / p
        N += 1

    # Analytic majorant for the full inclusive prefix T(q).
    # A(q) <= T(Q1) log(q)^(3/2)/sqrt(q) + K (q>=Q1).
    # K is the Appendix-A value; we recompute an outward upper enclosure.
    L0 = iv.log(iv.mpf(Q1))
    beta = iv.exp(-iv.euler) * (1 - 1/L0**2)
    F = (1 + iv.mpf('1.2762')/L0) * (1 + 1/L0)
    G = 1/(1 - 3/L0)
    factor = iv.sqrt((1+1/L0**2) * iv.mpf(Q1)/(Q1-1))
    K = factor * iv.exp(-iv.euler/2) * (1+iv.mpf('1.2762')/L0+F*G)
    L = iv.log(iv.mpf(X1))
    prefix = T_Q1 * iv.sqrt(L/iv.mpf(X1))
    if model == 'ellipse':
        constant = coeff
    elif design == 'disciplined':
        constant = 3 + Cf
    else:
        constant = iv.mpf(2)
    bound = coeff/iv.sqrt(beta) * (prefix + K/L) \
        + constant*iv.sqrt(L)/(iv.sqrt(beta)*iv.sqrt(X1)) \
        + (1 + iv.mpf('1.2762')/L)/L

    # All terms decay for q>=X1 (X1>=1e7 by default).
    analytic_ok = bound.b < iv.mpf(1).a
    finite_ok = last_fail is None
    print(f'{model} {design}: head last prime={head_prime}, window to {X1:,}')
    print(f'  primes checked: {N:,}; first after head: {first_after_head}')
    print(f'  max finite upper endpoint: {max_after[0]:.10f} at p={max_after[1]}')
    print(f'  last finite failure: {last_fail}')
    print(f'  Appendix-A K upper: {high(K):.9f}; beta lower: {low(beta):.9f}')
    print(f'  infinite-tail majorant upper at X1: {high(bound):.10f}')
    status = 'CERTIFIED' if finite_ok and analytic_ok else 'FAILED'
    print(f'  STATUS: {status} (rigorous outward-rounded interval checks)')
    return finite_ok and analytic_ok


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model', choices=['ellipse','triangle'])
    ap.add_argument('--design', choices=['disciplined','uniform8'], default='disciplined')
    ap.add_argument('--rho', default='0.6296271949192772')
    ap.add_argument('--p0', type=int, default=1000)
    ap.add_argument('--Q1', type=int, default=100000)
    ap.add_argument('--X1', type=int, default=10000000)
    ap.add_argument('--prec', type=int, default=128)
    ap.add_argument('--head-prime', type=int)
    a=ap.parse_args()
    raise SystemExit(0 if certify(a.model,a.design,a.rho,a.p0,a.Q1,a.X1,a.prec,a.head_prime) else 1)


if __name__ == '__main__':
    main()
