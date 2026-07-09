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
| `rect_exact.py` | The similar-rectangle model: guillotine scheduler, immature-stock measurement, and the two rigorous certificates behind Theorem 6.7 (`head`, `tail`). |
| `make_data.py` | Enumeration companion: produces the JSON data files consumed by the figure script (greedy triangle packing, integer dissection, adjacency chain, free-region decomposition). |
| `make_figures.py` | Renders the paper's figures from the JSON data (plus the two self-contained figures). |
| `head_certificate_sqrt2.json` | The 5133-prime interval certificate for the √2 head (Theorem 6.7), pre-generated. |

Generated on demand (not committed): `pieces_greedy.json`, `integers.json`, `pieces_chain_dfs.json`, `head_certificate.json`, and the figure `*.png` files.

---

## Requirements

Python 3.9+ and:

```
pip install numpy mpmath shapely matplotlib
```

* `mpmath` — 200-bit interval arithmetic for the rectangle certificates (`rect_exact.py`).
* `shapely` — exact polygon geometry for the triangle packer (`make_data.py`).
* `numpy`, `matplotlib` — sieving, statistics, and figures.

`rect_exact.py` needs only `numpy` + `mpmath`; the figure pipeline needs `shapely` + `matplotlib`.

---

## Quick start — reproduce the main theorem's certificates

```bash
# Finite head: place the first 5133 primes (q <= 5*10^4) in the sqrt2 rectangle,
# certify every placement fits in 200-bit interval arithmetic, and dump the full record.
python rect_exact.py head --q0 50000 --prec 200 --dump head_certificate_sqrt2.json

# Analytic tail: verify the immature-stock ratio rho(q) < 1 for all q > q0
# (finite windows + explicit Rosser-Schoenfeld/Dusart bound; Lemma 6.6 / Appendix A).
python rect_exact.py tail
```

Expected head summary: `CERTIFIED`, 5133 placements, reserve hit once (at *p*=3) with slack 0.092,
reserve aspect √2 to the interval, reserve+fringe area = R_{<q₀} exactly, minimum fit margin 2.4×10⁻³.
Expected tail summary: `CERTIFIED`, ρ(q) ≤ 0.923 < 1−1/q for every prime q > q₀ = 5·10⁴.

Together these are the two halves of the unconditional, computer-assisted proof of Theorem 6.7.

---

## `rect_exact.py` — the rectangle theorem

```
python rect_exact.py head     [--q0 50000] [--prec 200] [--reserve-frac 0.5] [--dump FILE]
python rect_exact.py tail     [--q0 50000] [--Q1 100000] [--X1 10000000]
python rect_exact.py ledger   [--nmax 700000]     # pure fattest-first guillotine run (float)
python rect_exact.py backstop [--nmax N]          # reserve-backstop scheduler run (float)
```

* **`head`** — runs the reserve-backstop scheduler for primes ≤ `q0` in `prec`-bit interval arithmetic, and
  certifies the handoff Theorem 6.5 requires: free region = a reserve rectangle similar to the container
  (area κ₀) + finitely many fringe rectangles, with every placement fitting by a certified positive margin.
  `--dump FILE` writes the full certificate (see schema below).
* **`tail`** — verifies ρ(q) = (3T(q)+r)/√(q·R_{<q}) < 1 for all q > q0 across three ranges:
  (q₀,10⁵] and (10⁵,10⁷] by direct evaluation, and (10⁷,∞) by the explicit analytic majorant of Appendix A
  (constants β, K derived from Rosser–Schoenfeld and Dusart). Prints the crossover (last q with ρ ≥ 1) and
  confirms it lies below q₀.

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
| `staircase_angle` | none | Staircase limit angle θ* = arctan(e^L) ≈ 56.29° (Theorem 3.1). |
| `rect_packing` | none | Exact guillotine packing of the √2 rectangle (Figure for §6). |
| `integers` | `integers.json` | Recursive least-prime-factor dissection; the sub-dissection depends only on a cell's largest prime factor, so e.g. [15] is an exact 1/√3-scaled copy of [5] (Proposition 4.2). |
| `vertex_stats` | `pieces_greedy.json` | Vertex statistics of the greedy triangle packing (side distances, angular coordinate, angle histogram). |

`make_data.py` also has `chain` (adjacency-chain DFS → `pieces_chain_dfs.json`) and `certificate`
(free-region decomposition → `head_certificate.json`); `make_figures.py` can render the corresponding
`packing`, `chain`, and `certificate` figures. These are **not** used by the current paper.

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
| `q0`, `prec_bits`, `n_primes` | head cutoff (5·10⁴), interval precision (200), primes placed (5133) |
| `reserve_min_slack` | min over reserve hits of √κ − √d (linear reserve safety), = 0.092 |
| `reserve_aspect` | `[lo, hi]` interval enclosing the reserve's aspect ratio (encloses √2) |
| `reserve_area` | `[lo, hi]` interval for the final reserve area κ₀ |
| `R_lt_q0` | `[lo, hi]` for R_{<q₀} = ∏_{p<q₀}(1−1/p) |
| `reserve_plus_fringe` | `[lo, hi]` for reserve + fringe areas (encloses `R_lt_q0` — area conservation) |
| `min_fit_margin` | minimum certified fit margin over all placements (≥ 0) |
| `reproduce` | the exact command that regenerates the file |
| `schema` | one-line description of a placement record |
| `placements` | array of per-prime records |

Each entry of `placements`:

```json
{
  "p": 3,
  "piece_sides": [[a_lo, a_hi], [b_lo, b_hi]],   // piece side lengths as 200-bit intervals (a >= b)
  "host_sides":  [[A_lo, A_hi], [B_lo, B_hi]],   // hosting rectangle side lengths as intervals
  "fit_margin":  0.077,                          // certified lower bound on min(A-a, B-b), >= 0
  "kind": "reserve"                              // "fringe" or "reserve" (a reserve hit)
}
```

Every `fit_margin ≥ 0` is a rigorous interval-arithmetic certificate that the piece fits its host.

---

## Result → command map

| Paper result | Command |
|--------------|---------|
| Theorem 6.7 (√2 packing), finite head | `python rect_exact.py head --q0 50000 --dump head_certificate_sqrt2.json` |
| Theorem 6.7, analytic tail (Lemma 6.6 / App. A) | `python rect_exact.py tail` |
| Theorem 3.1 (staircase angle) | `python make_figures.py staircase_angle` |
| Proposition 4.2 (self-similar dissection) | `python make_data.py integers && python make_figures.py integers` |
| §6 packing figure | `python make_figures.py rect_packing` |
| §7 vertex statistics | `python make_data.py greedy --nmax 355 && python make_figures.py vertex_stats` |

---

## Citing

If you use this code, please cite the paper. The certificate `head_certificate_sqrt2.json` and the finite-window
evaluations are the archival record behind the computer-assisted Theorem 6.7.

## License

MIT — see `LICENSE`.
