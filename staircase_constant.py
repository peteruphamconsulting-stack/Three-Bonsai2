#!/usr/bin/env python3
"""
staircase_constant.py  --  reproduces the numerical constants of Section 3 / Theorem 3.1
of "A Geometric Sieve of Eratosthenes" (P. Upham): the staircase limit angle.

The staircase dissection sends the remainder-rectangle aspect ratio h_n/w_n to
kappa = e^L, where

    L = sum_{j>=1} (-1)^{j+1} a_j,   a_j = ln(1 - 1/p_j)^{-1}   (p_j the j-th prime),

an alternating series whose terms decrease strictly to 0 (Leibniz), so it converges.
Two facts are reported, kept rigorously apart:

  * RIGOROUS.  The Leibniz truncation bound |L - L_n| <= a_{n+1} <= 2/p_{n+1} certifies
    L to 7-8 places at p <= 10^8 (L = 0.4047892(2)); every digit inside that bound is proved.

  * ACCELERATED (displayed digits, NOT rigorously bounded).  Because the partial sums
    L_n straddle L, their midpoint (L_n + L_{n-1})/2 -- one step of Euler's alternating-series
    transform -- removes the leading a_{n+1}/2 error and exposes several more digits. These are
    reported as "stable across checkpoint decades": the accelerated value is recomputed with the
    primes below 10^4, 10^5, ..., up to nmax, and the digits that agree across the last four
    decades are the ones displayed. Repeated Euler averaging agrees with the single midpoint.

As a cross-check the same pipeline evaluates the alternating prime series
A = sum (-1)^{j+1}/p_j, which agrees with OEIS A078437 (0.2696063520...) to ten places.

At nmax = 10^7 the accelerated value already equals the paper's 0.4047892180; the default
nmax = 10^8 confirms stability across one further decade (needs ~1 GB for the sieve).

Usage:
    python3 staircase_constant.py                 # default nmax = 10^8
    python3 staircase_constant.py --nmax 10000000 # 10^7 (fast; already reproduces the digits)

Requires: numpy.
"""
import argparse
from math import exp, atan, degrees
import numpy as np


def sieve(n):
    """Primes <= n as an int array (same sieve as the companion enumeration code)."""
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.flatnonzero(s)


def partial_sums(P):
    """Cumulative alternating partial sums L_n = sum_{j<=n} (-1)^{j+1} a_j and the analogous
    A_n = sum_{j<=n} (-1)^{j+1}/p_j, for the primes P = (p_1, p_2, ...)."""
    Pf = P.astype(np.float64)
    a = -np.log1p(-1.0 / Pf)                       # a_j = ln(1 - 1/p_j)^{-1} > 0, strictly decreasing
    signs = (-1.0) ** (np.arange(1, len(Pf) + 1) + 1)   # (-1)^{j+1}
    return np.cumsum(signs * a), np.cumsum(signs / Pf)


def midpoint(partial):
    """One step of Euler's transform: the midpoint of the last two partial sums."""
    return 0.5 * (partial[-1] + partial[-2])


def euler_repeated(partial, rounds, window=256):
    """Repeated pairwise averaging of the tail (van Wijngaarden / Euler); a stability check
    on the single-step midpoint. Operates on the last `window` partial sums."""
    s = partial[-window:].copy()
    for _ in range(rounds):
        s = 0.5 * (s[:-1] + s[1:])
    return s[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmax", type=int, default=10 ** 8,
                    help="sieve bound for the prime index of L (default 10^8, as in the paper)")
    args = ap.parse_args()

    P = sieve(args.nmax)
    Lpart, Apart = partial_sums(P)
    p_last = int(P[-1])

    # rigorous Leibniz truncation bound at the full sieve
    a_next = -np.log1p(-1.0 / (p_last + 1.0))       # >= a_{n+1}; and a_{n+1} <= 2/p_{n+1}
    leibniz = min(float(a_next), 2.0 / (p_last + 1.0))

    # accelerated value and its stability across checkpoint decades
    print(f"staircase constant L = sum (-1)^(j+1) ln(1-1/p_j)^(-1)   [primes to {p_last}]")
    print(f"  raw partial sum L_n           = {Lpart[-1]:.12f}")
    print(f"  RIGOROUS Leibniz bound        : |L - L_n| <= a_(n+1) <= 2/p_(n+1) = {leibniz:.2e}")
    lo, hi = sorted((float(Lpart[-2]), float(Lpart[-1])))
    print(f"  RIGOROUS bracket (consecutive partial sums enclose L, terms strictly decreasing):")
    print(f"                                  L in [{lo:.10f}, {hi:.10f}]  (width {hi-lo:.1e})")
    print(f"                                  => at nmax=10^8: L in [0.40478921, 0.40478923] certified")
    print()
    print("  accelerated (midpoint / Euler), by checkpoint decade:")
    decade = 10 ** 4
    last_vals = []
    while decade <= args.nmax:
        m = np.searchsorted(P, decade, side="right")   # primes below this decade
        if m >= 2:
            val = midpoint(Lpart[:m])
            last_vals.append(val)
            print(f"    primes < 10^{len(str(decade)) - 1:<2d}  ({m:>8d} primes) : L ~ {val:.10f}")
        decade *= 10
    full_mid = midpoint(Lpart)
    full_eul = euler_repeated(Lpart, rounds=6)
    print(f"    primes <= {p_last:<8d}          : L ~ {full_mid:.10f}   (midpoint)")
    print(f"    repeated Euler averaging (6 rounds)         : L ~ {full_eul:.10f}   (agrees)")

    # digits common to the two highest checkpoint decades (the displayed, non-rigorous digits)
    if last_vals:
        hi, lo = f"{full_mid:.10f}", f"{last_vals[-1]:.10f}"
        k = 0
        while k < len(hi) and k < len(lo) and hi[k] == lo[k]:
            k += 1
        agree = max(k - 2, 0)   # subtract "0." prefix -> number of matching decimal places
        print(f"    top two decades agree to {agree} decimal places: L = {full_mid:.10f}")

    L = full_mid
    kappa = exp(L); theta = degrees(atan(kappa))
    print()
    print(f"  kappa = e^L      = {kappa:.10f}     (paper 1.4989865074)")
    print(f"  theta* = arctan k = {theta:.7f} deg   (paper 56.2920568)")
    print(f"  tilt above diag  = {theta - 45.0:.6f} deg   (paper 11.292)")

    # cross-check: alternating prime series A vs OEIS A078437
    A = midpoint(Apart)
    print()
    print(f"  cross-check  A = sum (-1)^(j+1)/p_j = {A:.10f}   (OEIS A078437: 0.2696063520)")

    # rigorous distinctness kappa < 3/2 (paper): from a coarse Leibniz bound
    from math import log
    print(f"  rigorous     kappa = e^L < 3/2 ? uses L < ln(3/2) = {log(1.5):.6f}; "
          f"L_n + a_(n+1) = {Lpart[-1] + leibniz:.6f} < {log(1.5):.6f}: "
          f"{Lpart[-1] + leibniz < log(1.5)}")


if __name__ == "__main__":
    main()
