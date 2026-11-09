# Companion code — *A Geometric Sieve of Eratosthenes*

Reproduction code and certified data for the paper
**"A Geometric Sieve of Eratosthenes: Least-prime-factor densities, perfect packings, and the integers as a point set"** (P. Upham).

For each prime *p*, let *d_p* be the natural density of the integers whose least prime factor is *p*
(so *d_p* = (1/p)∏_{q<p}(1−1/q) and Σ_p *d_p* = 1). This code certifies the paper's three unconditional
perfect-packing results and generates its figures and empirical data:

* **Rectangle (§5).** The complete feasible interval for axis-parallel packings by similar rectangles is
  **1/√2 + 1/√6 ≤ r ≤ √6**. Sixteen interval-certified guillotine schedules of 1,282 pieces each cover
  the entire closed parameter interval; an arrangement-independent tail criterion places every later piece.
  Moreover every shifted prime state *s* ≥ 3 admits a perfect square packing, furnishing the recursive states.
* **Ellipse (§6).** The ellipse of aspect ratio ρ ≈ 0.6296 is perfectly packed by similar ellipses of areas
  *d_p* in two orientations — a 1650-piece certified head (exact rational Perram–Wertheim disjointness
  witnesses, interval containment) plus the **corrected orientation-averaged** eroded-volume capacity tail, and the ten further certified
  *state-head files* (nine required beyond state 2, plus the now-redundant state-31 certificate) that carry the recursive integer dissection into elliptical cells.
* **Triangle (§7).** The isosceles right triangle is perfectly packed by similar triangles of areas *d_p*
  in the eight 45° orientations — a 10906-piece certified head (exact-rational and 100-bit interval
  containment margins, separating-axis disjointness witnesses on a rational-edge frame) plus the same
  eroded-volume capacity tail, now with the triangle's mixed-area constant κ_T = 3/√2, and the certified
  family of **45 per-state heads** (states 2…197) that carry the recursive integer dissection into
  triangular cells. An exact 15-piece algebraic base, verified symbolically in ℚ(√2,…,√47), anchors the
  head's first pieces.
* **Staircase (§3).** The staircase limit-angle constant L and θ\* = arctan(e^L) ≈ 56.29°.

### Revised certification (November 8–9, 2026)

The rectangular proof now uses `rectangle_full_interval_certificate.py` instead of the
earlier reserve/immature-stock argument. Running the driver reproduced a **CERTIFIED** complete interval:
16 schedules, both exact algebraic endpoint contacts (2,563 strict remaining inequalities per endpoint),
maximum finite-tail upper value `0.9998320673924572` at prime `10487`, infinite majorant
`0.7616246752362302`, and the aligned-square state comparison for every prime state *s* ≥ 3.
It writes `rectangle_full_interval_certificate.json`.

The previous ellipse mixed-area coefficient `(rho + 1/rho)/2` was incorrect for
orthogonally oriented ellipses. The repaired proof uses
`kappa_star = 2 E(1-rho^4)/(pi*rho)` (complete elliptic integral of the second kind,
parameter convention) and its orientation average `kappa_bar = (1+kappa_star)/2`.
The analytic majorant now retains the complete prime-prefix contribution and
the `N_q/q` prime-count term. The finite window extends to `10,000,000` and is
checked by outward-rounded 128-bit interval arithmetic.

Two independent user-side runs of the patched scripts returned `STATUS: CERTIFIED`:

| Model | Head retained | Largest post-head criterion (upper endpoint) | Infinite-tail bound (upper endpoint) |
|---|---|---:|---:|
| Ellipse, two-orientation average | 1,650 pieces through 13,967 | 0.9769179523 at 14,011 | 0.5108274068 |
| Triangle, disciplined schedule | 10,906 pieces through 115,363 | 0.9999700506 at 115,399 | 0.6917952459 |

Each run tested `664,579` primes through `10^7` with **no finite post-head
failure**. The common explicit constants were `K <= 2.063589925` and
`beta >= 0.557223571`. The ellipse run used the exact decimal
`rho = 0.6296271949192772` from the existing head; it reported
`kappa_star ≈ 1.161361468300` and `kappa_bar ≈ 1.080680734150`.

The original geometric head JSON files and figure data are **unchanged**.
On November 9 the 1,650-piece ellipse container head was separately reverified on the author's Mac:
all `1,360,425` pairs passed (5,771 exact rational Perram–Wertheim witnesses),
the minimum exact PW slack was `0.000170927`, and the maximal containment support gap
was `−1.86368e−06`, strictly negative. The revised container-tail certificates are separate
from the existing geometric head checks. The current averaged ellipse state comparison makes
state 31 head-free; the older state-31 JSON remains valid but is no longer needed in the proof.

---

## Contents

### Rectangle proof (§5)

| File | Role |
|------|------|
| `rectangle_full_interval_certificate.py` | **Current main theorem driver:** certifies 16 parameter-interval schedules, both exact aspect-ratio endpoints, the mixed-orientation uniform tail, and every shifted square state; creates the summary JSON. |
| `rectangle_aspect_path_range.py` | Reconstructs fixed guillotine schedules and certifies all fit inequalities across intervals in *r*, including exact endpoint contacts. |
| `rectangle_aspect_sweep.py` | Prime enumeration and deterministic finite-head scheduling support. |
| `rectangle_aspect_recursive_certificate.py` | Verifies the aligned-square state heads and head-free state comparison, invoked by the main driver. |
| `rectangle_aspect_head_certificate.py` | Interval verifier used by the recursive-state driver. |
| `rectangle_full_interval_certificate.json` | Output of the current sharp-interval theorem driver (16 interval records, schedule hashes, and tail bounds). Regenerate with the main command below. |
| `rect_exact.py`, `head_certificate_sqrt2.json` | **Legacy √2-specific construction:** valid earlier guillotine/reserve experiment, retained for reproducibility and figures, but not the proof of the current sharp aspect-ratio theorem. |

### Ellipse proof (§6)

| File | Role |
|------|------|
| `ellipse_certify.py` | The geometric head certifier is unchanged in substance: exact rational Perram–Wertheim pair witnesses and interval containment. Its revised `tail` command calls `revised_tail_certificate.py`, uses the correct orientation-averaged mixed area, and certifies the post-head interval and infinite tail. The existing head ends at 13,967; this is a retained head cutoff, **not** the improved numerical criterion crossover. |
| `ellipse_build.py` | The head builder: beam search with spatial hashing, Sobol candidate fields, Nelder–Mead refinement gated on the exact margin the certifier will see, checkpoint/resume, and jam repair with escalating backup windows. `--state s` builds the state-*s* heads. Output JSON feeds `ellipse_certify.py` directly. |
| `ellipse_head_n1650.json` | The certified **container head** (state 2): 1650 pieces, primes 2…13967, exact dyadic centers, orientation pattern, ρ = 0.6296271949192772. |
| `head_s3.json` … `head_s31.json` | Ten additional existing certified **state-head files** (states 3, 5, 7, 11, 13, 17, 19, 23, 29, 31; sizes 1044, 755, 594, 474, 405, 321, 273, 209, 171, 115). Under the revised averaged criterion, only states 3 through 29 need these nine heads; **every state s ≥ 31 is head-free**. Thus ten heads total, including the container's state-2 head, suffice. The state-31 file may be retained as extra verification data. |
| `revised_tail_certificate.py` | Shared rigorous 128-bit interval certificate for the ellipse and triangle **container tails**, with a finite prime window through `10^7`, a full Appendix-A analytic majorant, and fail-closed exit status. |

### Triangle proof (§7)

| File | Role |
|------|------|
| `triangle_certify.py` | The rigorous triangle geometric head certifier is unchanged in substance. Its revised `tail --design disciplined` delegates to `revised_tail_certificate.py`, checking all post-head primes and the complete analytic majorant including `N_q/q`. The existing head ends at 115,363 and retains 10,906 pieces. |
| `triangle_build.py` | The head builder: corner-anchored greedy over 24 templates (8 orientations × 3 anchor vertices), contact-length (deepest-nestle) scoring, vectorized separating-axis overlap tests with a spatial grid, touch-then-nudge clearance snapped to the 2⁻⁴⁰ dyadic grid, escalating back-off jam repair, checkpoint/resume. `--state s` builds the state-*s* head (target auto-computed from the crossover); `--design disciplined` is the certified schedule. Output JSON feeds `triangle_certify.py` directly. |
| `run_triangle_states.py` | Drives the whole integer family: for each prime state 2…terminal it builds (via `triangle_build`) and certifies (via `triangle_certify`) one head, reuses any already-`CERTIFIED` head, and writes `states/tri_state_<s>.json` plus `states/manifest.json`. `--only`, `--keep-going`, `--certify-only` for subsets / re-verification. |
| `triangle_head_n10906.json` | The certified **container head** (state 2): 10906 pieces, primes 2…115363, disciplined schedule, rational-edge frame (`hypotenuse-on-x-axis`). Schema below. |
| `states/tri_state_<s>.json` | The **45 certified state heads** (states 2, 3, 5, …, 197; sizes 10906 down to 398, total 118593 pieces). Every state *s* ≥ 199 is covered by the tail from its first piece, so these 45 heads are exactly what the recursive dissection of ℕ into triangular cells requires. Same schema as the container head. |
| `states/manifest.json` | Family manifest: design, p0, per-state `{state, n, head, status}`, total piece count, terminal state, and `all_certified`. |
| `triangle_base_exact.py` | The exact algebraic base verifier (Remark "An exact algebraic base"): reconstructs the first 15 pieces (primes 2…47) in ℚ(√2,…,√47), proves each area = *d_p*, and decides all 135 containment sign-tests and 105 separating-edge certificates exactly — a stricter, fully symbolic check on the head's first pieces. |
| `_greedy.py` | The greedy placement core used by `triangle_base_exact.py` — extracted **verbatim** from `make_data.py` and frozen, so the exact verifier reproduces the same combinatorial choices without depending on the figure pipeline. The duplication with `make_data.py` is deliberate. |
| `base_combinatorics.json` | Frozen 15-piece base combinatorics (prime, anchor, rotation, per-edge direction, exact vertices); lets `triangle_base_exact.py --frozen` re-verify the base packing with no float greedy / GEOS dependence. |

### Staircase constant (§3)

| File | Role |
|------|------|
| `staircase_constant.py` | Reproduces the §3 numbers: the staircase constant L via midpoint (Euler) acceleration to 10⁸ with checkpoint decades and the rigorous Leibniz truncation bound, κ = e^L, θ\*, and the OEIS A078437 cross-check. Rigorous and accelerated digits are kept strictly apart. |

### Figures

| File | Role |
|------|------|
| `make_figures.py` | Renders `staircase_angle.png`, `rect_packing.png`, and `vertex_stats.png` (plus the legacy `integers`, `packing`, `chain`, `certificate` figures, no longer used by the paper). |
| `make_data.py` | Enumeration companion: produces `pieces_greedy.json` for the vertex-statistics figure (plus legacy `integers` / `chain` / `certificate` targets). |
| `rect_recursion_figure.py` | `rect_recursion.png` — the integers in the √2 rectangle, recursive least-prime-factor guillotine. Self-contained (the guillotine is exact and deterministic). |
| `ellipse_figure.py` | `ellipse_packing.png` — the certified 1650-piece ellipse head, coloured by scale, largest cells labelled by prime. Reads `ellipse_head_n1650.json`. |
| `ellipse_recursion_figure.py` | `ellipse_recursion.png` — the integers in the ellipse: the container head with the state-3/5/7 heads mapped inside the cells 3, 5, 7 by the exact similarity. Reads the head JSONs. |

Generated on demand (not committed): `pieces_greedy.json` and the figure `*.png` files.

---

## Requirements

Python 3.9+ and:

```
pip install numpy scipy mpmath sympy shapely matplotlib
```

Tested with (pin these for exact reproduction of the certificates and the frozen base record):
`numpy 2.4.4`, `scipy 1.x`, `sympy 1.14.0`, `mpmath 1.3.0`, `shapely 2.1.2` (GEOS 3.x), `matplotlib 3.x`, Python 3.11.

* `mpmath` — interval arithmetic for the rectangle certificates (200-bit), the ellipse containment
  certification, triangle containment margins (`triangle_certify.py`, 100-bit), and the revised
  ellipse/triangle tail inequalities (`revised_tail_certificate.py`, 128-bit).
* `sympy` — exact symbolic computation in ℚ(√2,…,√47) (`triangle_base_exact.py`) and exact prime/density
  bookkeeping for the ellipse and triangle scripts.
* `shapely` — floating-point polygon geometry (GEOS); used **only** to select packing combinatorics and to
  render figures. Every certified inequality is decided exactly (`sympy` / `Fraction`) or by
  directed-rounding interval arithmetic (`mpmath`), never by `shapely`.
* `scipy` — Sobol sampling and Nelder–Mead refinement in `ellipse_build.py` (search only; nothing certified
  depends on it).
* `numpy`, `matplotlib` — sieving, statistics, and figures.

Per script: `revised_tail_certificate.py` needs only `mpmath`; `rect_exact.py` needs only `numpy` + `mpmath`; `ellipse_certify.py` needs `numpy` + `mpmath` +
`sympy`; `ellipse_build.py` needs `numpy` + `scipy` + `sympy`; `triangle_certify.py` needs `numpy` + `mpmath`
+ `sympy`; `triangle_build.py` needs only `numpy`; `run_triangle_states.py` needs `numpy` + `sympy` (it imports
`triangle_build` and shells out to `triangle_certify`); `triangle_base_exact.py` needs `numpy` + `sympy` +
`mpmath` + `shapely` (and its `_greedy.py` core); the figure pipeline needs `shapely` + `matplotlib` +
`sympy`.

---

## Quick start — reproduce the certificates

### Rectangle (§5)

```bash
# Current sharp interval theorem: run from the folder containing the five
# rectangle_aspect*.py / rectangle_full_interval_certificate.py scripts.
python3 rectangle_full_interval_certificate.py
```

Expected concluding lines:

```text
CLOSED-INTERVAL COVER: CERTIFIED (16 schedules plus two exact contacts)
UNIFORM MIXED-ORIENTATION TAIL: CERTIFIED
ALL SHIFTED STATES s>=3: CERTIFIED (universal aligned-square family)
OVERALL: CERTIFIED COMPLETE ASPECT INTERVAL
Saved rectangle_full_interval_certificate.json
```

The generated JSON, the five scripts, and their exact version information belong in the archive.
For historical comparison, the earlier √2-only proof may still be reproduced with
`python3 rect_exact.py head --q0 50000 --prec 200 --dump head_certificate_sqrt2.json`
and `python3 rect_exact.py tail`. Those legacy commands are **not prerequisites** for the sharp theorem.

### Ellipse (§6)

```bash
# Certify the 1650-piece container head: exact sigma_2(q) areas, exact rational
# Perram-Wertheim disjointness for every close pair, interval containment for every piece.
python ellipse_certify.py certify --in ellipse_head_n1650.json

# Certify an existing shifted-state head (nine needed for states 3..29;
# the state-31 file is optional supplementary data):
python ellipse_certify.py certify --in head_s7.json

# Revised averaged tail: certified finite window to 10^7 and infinite majorant.
# The existing 1650-piece head through p=13967 is retained unchanged.
python ellipse_certify.py tail
```

Expected certify output: `STATUS: CERTIFIED` with all pairs disjoint (exact PW witnesses for close pairs,
rigorous distance pruning for far pairs) and a certified negative max support gap. The verified state-2 head
has `1,360,425` disjoint pairs, exact PW minimum slack `0.000170927`, and containment maximum support gap
`−1.86368e−06`. The container head
covers area Σ_{p≤13967} d_p = 1 − R_{<13997}; together with the `tail` crossover this is the
computer-assisted half of the ellipse theorem.
The certified post-head criterion has maximum upper endpoint `0.9769179523` at
prime `14011`; the analytic upper bound at `10^7` is `0.5108274068`.

### Triangle (§7)

```bash
# Certify the 10906-piece container head: exact sigma_2(q) areas, exact-rational or
# 100-bit interval containment margins, separating-axis disjointness witnesses.
python triangle_certify.py certify --in triangle_head_n10906.json --design disciplined

# Revised disciplined tail: 128-bit finite interval window to 10^7 and
# complete analytic majorant (including the prime-count term N_q/q).
python triangle_certify.py tail --design disciplined

# Build + certify the whole 45-state integer family in one pass (writes states/).
# Drop the container head in first so state 2 is reused, not rebuilt:
cp triangle_head_n10906.json states/tri_state_2.json
python run_triangle_states.py --design disciplined --p0 1000 --outdir states
```

Expected certify output: `STATUS: CERTIFIED` with all 98154 containment margins ≥ 0 (11 exact-zero leg
contacts), all pairs disjoint by separating-axis witnesses, least certified gap ≈ 3.5×10⁻¹⁸. Expected family
output: `45/45 states CERTIFIED, 118593 pieces`, `FAMILY CERTIFIED`, least gap over the family ≈ 4.2×10⁻¹⁸.
Together with the `tail` certificate, these supply the computational inputs for
the state-2 triangle packing theorem. The revised certified post-head maximum upper endpoint is
`0.9999700506` at prime `115399`; the analytic upper bound at `10^7` is
`0.6917952459`. The recursive integer corollary additionally uses the state-head
family and its analytic state comparisons.

```bash
# Optional stricter check: the exact algebraic base (first 15 primes) in Q(sqrt2..sqrt47).
python triangle_base_exact.py                       # float greedy fixes the combinatorics, then exact proof
python triangle_base_exact.py --frozen              # replay from base_combinatorics.json (no greedy/GEOS)
```

Expected: `CERTIFIED` with `135/135`, `105/105`, free area = R_{<53}; `--frozen` gives the same result from
the stored record alone.

### Staircase constant (§3)

```bash
python staircase_constant.py                        # nmax = 10^8 (paper's bound); --nmax 10000000 is fast
```

Expected: `L = 0.4047892180`, `kappa = 1.4989865074`, `theta* = 56.2920568°`, `A = 0.2696063520`
(OEIS A078437), with the Leibniz bracket certifying `L ∈ [0.40478921, 0.40478923]`.

---

## `rect_exact.py` — legacy √2-specific rectangle method

The following commands reproduce the earlier reserve/immature-stock construction,
retained for research provenance and the original figure. The current rectangle theorem
instead uses `rectangle_full_interval_certificate.py` as described in Quick start.

```
python rect_exact.py head     [--q0 50000] [--prec 200] [--reserve-frac 0.5] [--dump FILE]
python rect_exact.py tail     [--q0 50000] [--Q1 100000] [--X1 10000000] [--float]
python rect_exact.py ledger   [--nmax 700000]     # pure fattest-first guillotine run (float)
python rect_exact.py backstop [--nmax N]          # reserve-backstop scheduler run (float)
```

* **`head`** — runs the reserve-backstop scheduler for primes ≤ `q0` in `prec`-bit interval arithmetic, and
  certifies the handoff the buffer theorem requires: free region = a reserve rectangle similar to the
  container (area κ₀) + finitely many fringe rectangles, with every placement fitting by a certified
  positive margin. It also records the certified maturity-at-handoff facts (all 5135 inherited fringe
  rectangles contain the first post-head piece *q₁* = 50021). `--dump FILE` writes the full certificate
  (schema below).
* **`tail`** — verifies ρ(q) = (3T(q)+r)/√(q·R_{<q}) < 1 for all q > q0 across three ranges:
  (q₀,10⁵] and (10⁵,10⁷] by rigorous `mpmath` interval-arithmetic evaluation (each reported bound is a
  guaranteed upper bound; pass `--float` for the faster float64 path), and (10⁷,∞) by the explicit analytic
  majorant of Appendix A, whose constants and final assembly are themselves evaluated in interval
  arithmetic. Prints the crossover (last q with ρ ≥ 1) and confirms it lies below q₀.

The historical certificate compares side **lengths**; a positioned rendering of that
particular √2 packing is the `rect_packing` figure. Its reserve and maturity conditions are
not assumptions of the current arrangement-independent infinite continuation.

---

## `ellipse_certify.py` and `ellipse_build.py` — the ellipse theorem

```
python ellipse_certify.py certify [--in ellipse_head_n1650.json] [--cells 96]
python ellipse_certify.py tail    [--rho 0.6296271949192772]
python revised_tail_certificate.py ellipse [--X1 10000000] [--prec 128]
```

* **`certify`** — loads a head JSON, treating the stored binary64 centers as **exact dyadic rationals** and
  reconstructing every area share as the exact rational σ_s(q) = d_q/R_{<s} from *(state, n)* (their exact
  denominators are astronomical and are never parsed from the file). It then proves, fail-closed:
  * every close pair disjoint by an exact `Fraction` evaluation of the Perram–Wertheim contact function at
    a rational witness λ\* (far pairs pruned by a rigorous center-distance bound);
  * every piece contained by adaptive interval bisection (`mpmath.iv`) of the support gap in angle — a
    cell is accepted only when its interval upper bound is < 0.
* **`tail`** — checks the corrected averaged criterion
  `2*kappa_bar*(1+T_strict(q))/sqrt(q*R_<q)+N_q/q < 1`, where
  `T_strict(q)=sum_{p<q} sqrt(d_p)`,
  `kappa_star=2*E(1-rho^4)/(pi*rho)`, and `kappa_bar=(1+kappa_star)/2`.
  The elliptic integral is enclosed by a 240-term interval series. The finite
  prime window ends at `10^7` and the infinite tail uses the full, decreasing
  Appendix-A majorant. The previous claimed analytic evaluation at `3*10^4`
  is superseded. The original geometric head through 13,967 is retained.

The heads were discovered by `ellipse_build.py`:

```
python ellipse_build.py --target-n 1650 --out ellipse_head_n1650.json --checkpoint ckpt.json \
    [--resume HEAD_OR_CKPT.json] [--state 2] [--beam 6] [--sobol 3072] [--seed 7]
```

Beam search over Sobol candidate fields with a spatial hash (each candidate tested against O(1)
neighbours), Nelder–Mead refinement, and acceptance gated on the **exact margin the certifier will see**;
checkpoints every 25 pieces; on a jam it reshuffles the field, then backs up and re-places a window of
recent pieces with escalating window sizes. `--state s` builds the state-*s* heads (relative areas
σ_s(q)). The builder is a numerical **discovery** tool; all rigor lives in `ellipse_certify.py`.

### Ellipse head JSON schema

| Key | Meaning |
|-----|---------|
| `format`, `status` | `"ellipse-binary-head"`; status is advisory — run the certifier |
| `state`, `n`, `q_last` | least admissible child prime *s* (2 = container), piece count, last prime placed |
| `verification.aspect_ratio_b_over_a` | ρ (exact rational when parsed as a decimal literal) |
| `verification.orientation_pattern_degrees` | 0/90 per piece |
| `verification.centers` | actual centers (x, y) — exact dyadic rationals |
| `verification.centers_normalized` | (x, u) with y = ρu — the builder's working frame |
| `verification.area_shares_float` | float echoes of σ_s(q), matched (not trusted) by the certifier |
| `verification.sum_head_area_float` | float head area |

---

## `triangle_certify.py`, `triangle_build.py`, `run_triangle_states.py` — the triangle theorem

```
python triangle_certify.py certify [--in triangle_head_n10906.json] [--design disciplined] [--prec 100] [--state 2]
python triangle_certify.py tail    [--design disciplined] [--p0 1000] [--X 10000000] [--state 2]
python revised_tail_certificate.py triangle [--X1 10000000] [--prec 128]
python triangle_build.py [--state 2] [--design disciplined] [--p0 1000] [--target-n 0] [--out FILE] \
    [--checkpoint FILE] [--ckpt-every 500] [--resume FILE] [--seed 7] [--time-cap SECS]
python run_triangle_states.py [--design disciplined] [--p0 1000] [--outdir states] \
    [--only 197,193,191] [--max-state 400] [--keep-going] [--certify-only]
```

* **`certify`** — loads a head JSON, treats the stored dyadic anchor coordinates as **exact rationals**, and
  reconstructs every area share as the exact σ_s(q) = d_q/R_{<s} from *(state, n)* (never parsed from the
  file). It then proves, fail-closed:
  * every piece contained by a per-vertex edge margin against the three rational frame edges — decided by
    **exact `Fraction` arithmetic** where the vertex's radical coefficient vanishes (the intended leg /
    hypotenuse contacts, which certify as exact zeros), and by directed-rounding `mpmath.iv` interval
    arithmetic otherwise;
  * every pair disjoint by a separating-axis witness among the eight directions *k*·45° — these contain
    every edge normal of every piece — a certified positive gap being the witness (far pairs pruned by a
    padded bounding box), in exact parallel to the ellipse's rational Perram–Wertheim witnesses.
* **`tail`** — the corrected disciplined eroded-volume criterion with the mixed-area
  table *V_p* and κ_T = 3/√2. It retains the existing head through 115,363,
  rigorously checks every later prime through `10^7`, and uses the complete
  analytic bound **including** `N_q/q`. The revised command certifies **state 2**;
  `--state s` for `s>2` is not presently supported by the new tail verifier and
  fails explicitly rather than displaying an unsupported certificate.
* **`build`** — the numerical **discovery** tool: corner-anchored greedy over the 24 templates, deepest-nestle
  scoring, separating-axis overlap with a spatial grid, touch-then-nudge clearance snapped to the 2⁻⁴⁰ grid,
  escalating jam repair, checkpoint/resume. State 2 is seeded with the exact canonical half-split (*p* = 2);
  all rigor lives in `triangle_certify.py`.
* **`run_triangle_states`** — builds and certifies the full family, reusing any head already marked
  `CERTIFIED` (drop the container head at `states/tri_state_2.json` to reuse it). Writes one head per state and
  `states/manifest.json`.

### Triangle head JSON schema

| Key | Meaning |
|-----|---------|
| `format` | `"triangle-head"` |
| `frame` | `"hypotenuse-on-x-axis"` — the rational-edge frame (base triangle (0,0),(1,0),(0,1); edges *y* ≥ 0, *x*+*y* ≤ 1, *y*−*x* ≤ 1) |
| `design`, `p0` | schedule (`disciplined`) and its free-orientation cutoff (1000) |
| `state`, `n`, `q_last` | least admissible child prime *s* (2 = container), piece count, last prime placed |
| `note` | advisory (`"BUILT (verify with triangle_certify.py certify)"`); status is **not** asserted — run the certifier |
| `pieces` | array of *n* records `[p, k, a, ax, ay]` |

Each piece record is `[p, k, a, ax, ay]`: prime *p*; orientation *k* ∈ {0,…,7} (rotation by *k*·45°);
anchor-vertex index *a* ∈ {0,1,2}; and the dyadic anchor point (*ax*, *ay*). The piece's vertices are
`anchor + L_p·(Rot_k(base) − Rot_k(base)[a])` with L_p = √(2 *d_p*) and base = (0,0),(1,0),(0,1) — so the whole
head is reconstructed exactly from these five numbers per piece, and the certifier recomputes areas and
contacts symbolically rather than trusting any stored geometry.

### `states/manifest.json` schema

| Key | Meaning |
|-----|---------|
| `design`, `p0` | schedule and cutoff for the family |
| `states` | array of `{ "state": s, "n": pieces, "head": "tri_state_s.json", "status": "CERTIFIED" }` |
| `total_pieces` | sum over the family (118593 for the disciplined *p0* = 1000 family) |
| `terminal_state` | first state with an empty head (199 — tail-covered from its first piece) |
| `all_certified` | `true` only if every listed state certified |

---

## Figures

The paper uses **six** figures. Three come from `make_figures.py`; three from the dedicated scripts.

```bash
# data for the vertex-statistics figure (a few minutes; the greedy is the slow step)
python make_data.py greedy --nmax 355                 # -> pieces_greedy.json

# the six paper figures
python make_figures.py staircase_angle rect_packing vertex_stats
python rect_recursion_figure.py    --out rect_recursion.png
python ellipse_figure.py           --in ellipse_head_n1650.json --out ellipse_packing.png
python ellipse_recursion_figure.py --container ellipse_head_n1650.json \
    --state3 head_s3.json --state5 head_s5.json --state7 head_s7.json --out ellipse_recursion.png
```

| Figure | Data needed | Content |
|--------|-------------|---------|
| `staircase_angle.png` | none | Staircase limit angle θ\* = arctan(e^L) ≈ 56.29° (§3). For the authoritative §3 numbers and their rigorous/accelerated split, run `staircase_constant.py` (the figure uses the same midpoint-accelerated *L*). |
| `rect_recursion.png` | none | The integers in the √2 rectangle: recursive least-prime-factor guillotine, exact self-similarity (§4). |
| `ellipse_recursion.png` | 4 head JSONs | The integers in the ellipse: the state-3/5/7 heads mapped into the cells 3, 5, 7 by the exact similarity (§6, Part B). |
| `rect_packing.png` | none | Exact guillotine packing of the √2 rectangle (§5). |
| `ellipse_packing.png` | `ellipse_head_n1650.json` | The certified 1650-piece ellipse head, coloured by scale (§6). |
| `vertex_stats.png` | `pieces_greedy.json` | Vertex statistics of the 355-piece greedy triangle packing (§8). |

`make_data.py` also has legacy `integers`, `chain`, and `certificate` targets, and `make_figures.py` the
corresponding `integers`, `packing`, `chain`, and `certificate` figures. These are **not** used by the
current paper (the old `integers` figure is superseded by `rect_recursion_figure.py`). Note that
`make_data.py certificate` is a *floating-point* free-region decomposition (visualization only); the exact,
rigorous certification of the triangle base packing is `triangle_base_exact.py`, not this command.

**Determinism note.** The triangle packer is a faithful re-implementation of the corner-anchored greedy
described in the paper; the exact placement sequence and piece count depend on the placement rule and
floating tolerances, so a run may differ slightly from the figures in the manuscript. The geometry, the
self-similarity, and the qualitative statistics are stable. Knobs are exposed in
`make_data.py::integer_dissection` (`min_rel`, `min_area`, `max_level`, `level1`) and `greedy_pack`
(`time_cap`).

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
| `handoff` | object recording the certified mature handoff (see sub-table) |
| `reproduce` | the exact command that regenerates the file |
| `schema` | one-line description of a placement record + the handoff block |
| `placements` | array of per-prime records (5132 entries: the primes 3–49999; piece 2 is the implicit half-split) |

The `handoff` object (every inherited fringe rectangle already contains the first post-head piece *q₁*):

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
| Staircase limit angle (§3, Theorem `thm:angle`) | `python staircase_constant.py` |
| Self-similar dissection figures (§4, Prop. `prop:markov`) | `python rect_recursion_figure.py` and `python ellipse_recursion_figure.py …` |
| Sharp rectangle aspect-ratio classification (§5, Theorem `thm:rect`), sixteen heads and uniform tail | `python3 rectangle_full_interval_certificate.py` |
| Shifted square-state continuation (§5) | Invoked automatically by `rectangle_full_interval_certificate.py` (or `python3 rectangle_aspect_recursive_certificate.py --r 1 --square-states`) |
| Earlier √2-only reserve construction (legacy) | `python3 rect_exact.py head --q0 50000 --dump head_certificate_sqrt2.json` and `python3 rect_exact.py tail` |
| Ellipse certified head (§6, Theorem `thm:ellhead`) | `python ellipse_certify.py certify --in ellipse_head_n1650.json` |
| Ellipse tail crossover (Lemma `lem:ellthresh`, Theorem `thm:ellipse`) | `python ellipse_certify.py tail` |
| Required ellipse state heads (§6, Prop. `prop:states`) | `python ellipse_certify.py certify --in head_sN.json` for N ∈ {3,5,7,11,13,17,19,23,29}; state 31 is head-free |
| Integers in elliptical cells (§6, Theorem `thm:ellintegers`) | state-2 head plus nine required shifted-state heads; `python ellipse_recursion_figure.py …` for illustration |
| Triangle certified head (§7, Theorem `thm:trihead`) | `python triangle_certify.py certify --in triangle_head_n10906.json --design disciplined` |
| Triangle tail crossover (Prop. `prop:tricross`, Theorem `thm:tripack`) | `python triangle_certify.py tail --design disciplined` |
| 45 state heads (§7, Prop. `prop:tristates`) | `python run_triangle_states.py --design disciplined --outdir states` |
| Integers in triangular cells (§7, Theorem `thm:triintegers`) | the 45 certified heads (`states/manifest.json` → `all_certified: true`) |
| Triangle exact algebraic base (§7, Remark `rem:tribase`) | `python triangle_base_exact.py` (or `--frozen`) |
| §5 packing figure | `python make_figures.py rect_packing` |
| §6 packing figure | `python ellipse_figure.py --in ellipse_head_n1650.json` |
| §8 vertex statistics | `python make_data.py greedy --nmax 355 && python make_figures.py vertex_stats` |

---

## Citing

If you use this code, please cite the paper. This paper is archived at Zenodo:
[https://doi.org/10.5281/zenodo.21283480](https://doi.org/10.5281/zenodo.21283480) (DOI `10.5281/zenodo.21283480`).

The principal rectangular archival input is `rectangle_full_interval_certificate.json`
together with its five Python verification scripts. The earlier `head_certificate_sqrt2.json`
is retained as legacy documentation. The ellipse archive contains `ellipse_head_n1650.json`
and `head_s3.json`–`head_s31.json` (with `head_s31.json` now optional for the theorem).
The triangle archive contains `triangle_head_n10906.json` and 45 state heads under `states/`
with their `manifest.json`. These data, the revised tail scripts, and the finite-window
interval evaluations support the three computer-assisted packing theorems.

## License

MIT — see `LICENSE`.
