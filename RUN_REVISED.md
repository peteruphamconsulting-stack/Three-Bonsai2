# Re-certified packing tails — November 2026

Place these three files together in the repository root, replacing the existing
`ellipse_certify.py` and `triangle_certify.py` while **adding**
`revised_tail_certificate.py`:

- `ellipse_certify.py`
- `triangle_certify.py`
- `revised_tail_certificate.py`

The geometry/head builders, existing head JSON certificates and figure scripts
do **not** need regeneration for the revised ellipse orientation-averaged proof.
The triangle head geometry is likewise unchanged.

## Run on a Mac (single commands, no additional flags)

```bash
python3 ellipse_certify.py tail
python3 triangle_certify.py tail --design disciplined
```

The two commands verify every prime beyond the original head cutoff through
`10,000,000` using 128-bit outward-rounded `mpmath.iv` arithmetic and verify
the infinite remainder by the complete Appendix-A majorant, including
`N_q/q`. The run may take several minutes. The status is `CERTIFIED` only if
**every** finite prime passes and the analytic bound is strictly below one.
If a criterion fails, the command exits nonzero and prints the last failure;
this does not demonstrate a packing impossibility.

The default ellipse aspect is the exact decimal used by the existing head
builder (`0.6296271949192772`). The source manuscript rounds this to
`0.629627`; if the precise head record uses a different aspect, use
`--rho` with that exact decimal on the standalone `revised_tail_certificate.py`.

### Independent invocation

```bash
python3 revised_tail_certificate.py ellipse
python3 revised_tail_certificate.py triangle
```

The original head verifiers remain available and unchanged in substance:

```bash
python3 ellipse_certify.py certify --in ellipse_head_n1650.json
python3 triangle_certify.py certify --in triangle_head_n10906.json --design disciplined
```

## What changed

1. Corrected orthogonal ellipse mixed area using a rigorous 240-term elliptic-integral series enclosure.
2. Replaced the false worst-case `kappa=(rho+1/rho)/2` with the valid orientation average `(1+kappa_star)/2`.
3. Used the **strict** prime prefix `sum_{p<q} sqrt(d_p)` for finite placeability.
4. Replaced the unproved `A(q)<=A(100000)` assumption with the complete Appendix-A upper bound.
5. Included the `N_q/q` term in triangle and ellipse analytic bounds.
6. Added fail-closed status and nonzero process exits on insufficient certificates.

## Outstanding issues

The rectangle's finite `head` code checks geometric fit and the mature-fringe
handoff. It does not prove the *infinite* immature-stock charging lemma. That
requires a separate analytic invariant for repeated subdivision and reserve
lanes; no amount of head regeneration alone closes it. `rect_exact.py` is
therefore **not replaced** by this patch.

The revised tail checker currently proves the container-state (`s=2`) ellipse
and triangle tails. Original state-head JSON files remain relevant. A complete
independent audit of all state-specific analytic inequalities and of the
rectangle descendant-charging invariant remains a manuscript-level task.

The original repository README's historical crossover discussion should be
updated after the new certificate outputs are obtained.
