# Companion code — *A Geometric Sieve of Eratosthenes*

Reproduction code and certified data for the paper
**"A Geometric Sieve of Eratosthenes: Least-prime-factor densities, perfect packings, and the integers as a point set"** (P. Upham).

For each prime *p*, let *d_p* be the natural density of the integers whose least prime factor is *p*
(so *d_p* = (1/p)∏_{q<p}(1−1/q) and Σ_p *d_p* = 1). This code (i) **certifies** the paper's unconditional
result — that the aspect-ratio-√2 rectangle is perfectly packed by similar rectangles of areas *d_p* — and
(ii) **generates** the figures and the empirical data for the triangle model.

---

## Contents

| File | Role |
|------|------|
| `rect_exact.py` | The similar-rectangle model: guillotine scheduler, immature-stock measurement, and the two rigorous certificates behind Theorem 5.7 (`head`, `tail`; the `tail` finite windows and analytic-tail constants now run in `mpmath` interval arithmetic, and `head` derives its `CERTIFIED` status from the recorded invariants rather than asserting it). |
| `triangle_base_exact.py` | The exact Proposition 6.3 verifier (with its `_greedy.py` core): reconstructs the 15-piece base packing in ℚ(√2,…,√47), proves each area = *d_p*, decides all 135 containment sign-tests and 105 separating-edge certificates exactly, and certifies that the free region contains a disk of radius ≥ 0.106 — its centre proven (in 200-bit interval arithmetic) to lie in the free region and its radius a floored interval lower bound on the centre's distance to the boundary — so the inradius is ≥ 0.106. |
| `make_data.py` | Enumeration companion: produces the JSON data files consumed by the figure script (greedy triangle packing, integer dissection, adjacency chain, free-region decomposition). |
| `make_figures.py` | Renders the paper's figures from the JSON data (plus the two self-contained figures). |
| `staircase_constant.py` | Reproduces the Section 3 / Theorem 3.1 numbers: the staircase constant *L* via midpoint (Euler) acceleration to 10⁸ with checkpoint decades and the rigorous Leibniz truncation bound, κ = e^L, θ*, and the OEIS A078437 cross-check. |
| `base_combinatorics.json` | Frozen 15-piece base combinatorics (prime, anchor, rotation, per-edge direction, and exact vertices); lets `triangle_base_exact.py --frozen` re-verify Proposition 6.3 with no float greedy / GEOS dependence. |
| `head_certificate_sqrt2.json` | The 5133-prime interval certificate for the √2 head (Theorem 5.7), pre-generated. |

Generated on demand (not committed): `pieces_greedy.json`, `integers.json`, `pieces_chain_dfs.json`, `head_certificate.json`, and the figure `*.png` files.

---

## Requirements

Python 3.9+ and:

```
pip install numpy mpmath sympy shapely matplotlib
```

Tested with (pin these for exact reproduction of the certificates and the frozen base record):
`numpy 2.4.4`, `sympy 1.14.0`, `mpmath 1.3.0`, `shapely 2.1.2` (GEOS 3.x), `matplotlib 3.x`, Python 3.11.

* `mpmath` — 200-bit interval arithmetic for the rectangle certificates (`rect_exact.py`) and the
  inscribed-disk certification in `triangle_base_exact.py`.
* `sympy` — exact symbolic computation in ℚ(√2,…,√47) for the Proposition 6.3 verifier (`triangle_base_exact.py`).
* `shapely` — floating-point polygon geometry (GEOS); used only to select the packing combinatorics and to render figures. Every certified inequality is decided exactly in ℚ(√2,…,√47) (`sympy`) or by directed-rounding interval arithmetic (`mpmath`), never by `shapely`.
* `numpy`, `matplotlib` — sieving, statistics, and figures.

`rect_exact.py` needs only `numpy` + `mpmath`; `triangle_base_exact.py` needs `numpy` + `sympy` + `mpmath`
+ `shapely` (and its `_greedy.py` core module); the figure pipeline needs `shapely` + `matplotlib`.

---

## Quick start — reproduce the main theorem's certificates

```bash
# Finite head: place the first 5133 primes (q <= 5*10^4) in the sqrt2 rectangle,
# certify every placement fits in 200-bit interval arithmetic, and dump the full record.
python rect_exact.py head --q0 50000 --prec 200 --dump head_certificate_sqrt2.json

# Analytic tail: verify the immature-stock ratio rho(q) < 1 for all q > q0
# (finite windows in interval arithmetic + explicit Rosser-Schoenfeld/Dusart bound; Lemma 5.6 / Appendix A).
python rect_exact.py tail

# Exact base packing (Proposition 6.3): 15 areas = d_p, 135 containment + 105 disjointness certificates
# in Q(sqrt2,...,sqrt47), free area = R_<50, and a certified inscribed disk of radius >= 0.106.
python triangle_base_exact.py                       # float greedy fixes the combinatorics, then exact/interval proof
python triangle_base_exact.py --frozen              # replay from base_combinatorics.json (no greedy/GEOS)
python triangle_base_exact.py --emit-record base_combinatorics.json   # (re)freeze the combinatorics

# Staircase limit angle (Theorem 3.1): L, kappa = e^L, theta*, the OEIS A078437 cross-check,
# and the rigorous Leibniz truncation bound, by midpoint (Euler) acceleration with checkpoint decades.
python staircase_constant.py                        # nmax = 10^8 (paper's bound); --nmax 10000000 is fast
```

Expected: `triangle_base_exact.py` prints `CERTIFIED` with `135/135`, `105/105`, free area = R_{<50},
and inscribed-disk radius ≥ 0.106 (centre certified in the free region); `--frozen` gives the same result
from the stored record alone. `staircase_constant.py` prints `L = 0.4047892180`, `kappa = 1.4989865074`,
`theta* = 56.2920568°`, and `A = 0.2696063520` (OEIS A078437), with the Leibniz bound proving `L = 0.4047892(2)`.

Expected head summary: `CERTIFIED`, `placements certified: 5133` (all primes 2…49999; piece 2 is the
implicit left-half split), reserve hit once (at *p*=3) with slack 0.092, 5135 fringe rectangles, reserve
aspect √2 to the interval, reserve+fringe area = R_{<q₀} exactly, minimum fit margin 2.4×10⁻³. The dumped
JSON records the 5132 explicit placements (primes 3–49999); `n_primes = 5133` counts piece 2 as well.
Expected tail summary: `CERTIFIED`, ρ(q) ≤ 0.923 < 1−1/q for every prime q > q₀ = 5·10⁴.

Together these are the two halves of the unconditional, computer-assisted proof of Theorem 5.7.

---

## `rect_exact.py` — the rectangle theorem

```
python rect_exact.py head     [--q0 50000] [--prec 200] [--reserve-frac 0.5] [--dump FILE]
python rect_exact.py tail     [--q0 50000] [--Q1 100000] [--X1 10000000] [--float]
python rect_exact.py ledger   [--nmax 700000]     # pure fattest-first guillotine run (float)
python rect_exact.py backstop [--nmax N]          # reserve-backstop scheduler run (float)
```

* **`head`** — runs the reserve-backstop scheduler for primes ≤ `q0` in `prec`-bit interval arithmetic, and
  certifies the handoff Theorem 5.5 requires: free region = a reserve rectangle similar to the container
  (area κ₀) + finitely many fringe rectangles, with every placement fitting by a certified positive margin.
  It also records the certified maturity-at-handoff facts of Lemma 5.6A (the 5135 inherited fringe rectangles
  all contain the first post-head piece *q₁* = 50021).
  `--dump FILE` writes the full certificate (see schema below).
* **`tail`** — verifies ρ(q) = (3T(q)+r)/√(q·R_{<q}) < 1 for all q > q0 across three ranges:
  (q₀,10⁵] and (10⁵,10⁷] by rigorous `mpmath` interval-arithmetic evaluation (each reported bound is a
  guaranteed upper bound; pass `--float` to fall back to the faster float64 path), and (10⁷,∞) by the explicit
  analytic majorant of Appendix A whose constants (β, K, T(Q₁)) and final assembly are themselves evaluated in
  interval arithmetic (Euler's γ enclosed rigorously; Dusart and Mertens-factor constants entered as outward-rounded
  intervals), so the reported tail bound ρ < 0.614 is a guaranteed upper bound.
  Prints the crossover (last q with ρ ≥ 1) and confirms it lies below q₀.

The certificate is a comparison of side **lengths**, so it is position-independent; a positioned rendering of
the same packing is the `rect_packing` figure.

---

## `make_data.py` and `make_figures.py` — the figures

The paper uses **four** figures. Two are self-contained; two need data.

```bash
# 1) data (a few minutes; the greedy triangle packing is the slow step)
python make_data.py greedy --nmax 355     # -> pieces_greedy.json
python make_data.py integers              # -> integers.json

# 2) figures
python make_figures.py staircase_angle integers rect_packing vertex_stats
```

| Figure | Data needed | Content |
|--------|-------------|---------|
| `staircase_angle` | none | Staircase limit angle θ* = arctan(e^L) ≈ 56.29° (Theorem 3.1). For the authoritative §3 numbers and their rigorous/accelerated split, run `staircase_constant.py` (the figure uses the same midpoint-accelerated *L*). |
| `rect_packing` | none | Exact guillotine packing of the √2 rectangle (Figure for §5). |
| `integers` | `integers.json` | Recursive least-prime-factor dissection; the sub-dissection depends only on a cell's largest prime factor, so e.g. [15] is an exact 1/√3-scaled copy of [5] (Proposition 4.2). |
| `vertex_stats` | `pieces_greedy.json` | Vertex statistics of the greedy triangle packing (side distances, angular coordinate, angle histogram). |

`make_data.py` also has `chain` (adjacency-chain DFS → `pieces_chain_dfs.json`) and `certificate`
(free-region decomposition → `head_certificate.json`); `make_figures.py` can render the corresponding
`packing`, `chain`, and `certificate` figures. These are **not** used by the current paper. Note that
`make_data.py certificate` is a *floating-point* free-region decomposition (visualization only); the exact,
rigorous certification of Proposition 6.3 is `triangle_base_exact.py` (see below), not this command.

**Determinism note.** The triangle packer is a faithful re-implementation of the corner-anchored greedy
described in the paper; the exact placement sequence and piece count depend on the placement rule and floating
tolerances, so a run may differ slightly from the figures in the manuscript. The geometry, the self-similarity,
and the qualitative statistics are stable. Knobs are exposed in
`make_data.py::integer_dissection` (`min_rel`, `min_area`, `max_level`, `level1`) and
`greedy_pack` (`time_cap`).

---

## `head_certificate_sqrt2.json` — schema

Top level:

| Key | Meaning |
|-----|---------|
| `model` | `"sqrt2 rectangle, r=sqrt(2)"` |
| `status` | `"CERTIFIED"` **only if every flag in `invariants` holds** — the status is *derived* from the invariants (fail-closed), never assigned |
| `q0`, `prec_bits`, `n_primes` | head cutoff (5·10⁴), interval precision (200), primes placed (5133, counting piece 2) |
| `invariants` | object of the derived gate flags: `placements_fit`, `reserve_aspect_encloses_sqrt2`, `reserve_plus_fringe_encloses_R`, `all_fringe_mature`, `both_excesses_positive`, `prime_count_ok`, `handoff_prime_ok` |
| `fringe_count` | number of fringe rectangles held at handoff (5135) |
| `reserve_min_slack` | min over reserve hits of √κ − √d (linear reserve safety), = 0.092 |
| `reserve_aspect` | decimal-string pair `[lo, hi]` enclosing the reserve's aspect ratio (encloses √2) |
| `reserve_area` | decimal-string pair `[lo, hi]` for the final reserve area κ₀ |
| `R_lt_q0` | decimal-string pair `[lo, hi]` for R_{<q₀} = ∏_{p<q₀}(1−1/p) |
| `reserve_plus_fringe` | decimal-string pair `[lo, hi]` for reserve + fringe areas (encloses `R_lt_q0` — area conservation) |
| `min_fit_margin` | minimum certified fit margin over all placements, as a single floored decimal string (≥ 0) |
| `handoff` | object recording the certified mature handoff of Lemma 5.6A (see sub-table) |
| `reproduce` | the exact command that regenerates the file |
| `schema` | one-line description of a placement record + the handoff block |
| `placements` | array of per-prime records (5132 entries: the primes 3–49999; piece 2 is the implicit half-split) |

The `handoff` object (Lemma 5.6A — every inherited fringe rectangle already contains the first post-head piece *q₁*):

| `handoff` key | Meaning |
|-----|---------|
| `next_prime` | the first prime past the head, *q₁* = 50021 |
| `next_piece_sides` | `{ "long_a": [lo, hi], "short_b": [lo, hi] }` — the *q₁*-piece side lengths (decimal-string endpoints) |
| `fringe_count` | fringe rectangles inherited at handoff (5135; echoes the top-level key) |
| `all_fringe_mature` | `true` — every inherited fringe rectangle contains the *q₁*-piece |
| `minimum_long_side_excess` | floored decimal string — min over fringe of (long side − *a_{q₁}*), ≥ 1.19758×10⁻³ |
| `minimum_short_side_excess` | floored decimal string — min over fringe of (short side − *b_{q₁}*), ≥ 1.76847×10⁻⁶ |

Interval endpoints throughout the file are serialized as **outward-rounded decimal-string pairs `[lo, hi]`**
(not a single collapsed float), so each interval keeps its certified width; one-sided certified lower bounds
(fit margins, side-excesses) are serialized as a single floored decimal string.

Each entry of `placements`:

```json
{
  "p": 3,
  "piece_sides": [["a_lo", "a_hi"], ["b_lo", "b_hi"]],   // decimal-string interval endpoints (a >= b)
  "host_sides":  [["A_lo", "A_hi"], ["B_lo", "B_hi"]],   // hosting rectangle side-length intervals
  "fit_margin":  "0.07715368...",                        // floored decimal string: certified lower bound on min(A-a, B-b), >= 0
  "kind": "reserve"                                       // "fringe" or "reserve" (a reserve hit)
}
```

Every `fit_margin ≥ 0` is a rigorous interval-arithmetic certificate that the piece fits its host.

---

## Result → command map

| Paper result | Command |
|--------------|---------|
| Theorem 5.7 (√2 packing), finite head | `python rect_exact.py head --q0 50000 --dump head_certificate_sqrt2.json` |
| Theorem 5.7, analytic tail (Lemma 5.6 / App. A) | `python rect_exact.py tail` |
| Theorem 3.1 (staircase angle) | `python make_figures.py staircase_angle` |
| Proposition 4.2 (self-similar dissection) | `python make_data.py integers && python make_figures.py integers` |
| Proposition 6.3 (exact base packing) | `python triangle_base_exact.py` |
| §5 packing figure | `python make_figures.py rect_packing` |
| §7 vertex statistics | `python make_data.py greedy --nmax 355 && python make_figures.py vertex_stats` |

---

## Citing

If you use this code, please cite the paper. This repository is archived at Zenodo:
[https://doi.org/10.5281/zenodo.21283480](https://doi.org/10.5281/zenodo.21283480) (DOI `10.5281/zenodo.21283480`).

The certificate `head_certificate_sqrt2.json` and the finite-window evaluations are the archival record
behind the computer-assisted Theorem 5.7.

## License

MIT — see `LICENSE`.
