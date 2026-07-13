 # Companion code — *A Geometric Sieve of Eratosthenes*

Reproduction code and certified data for the paper
**"A Geometric Sieve of Eratosthenes: Least-prime-factor densities, perfect packings, and the integers as a point set"** (P. Upham).

For each prime *p*, let *d_p* be the natural density of the integers whose least prime factor is *p*
(so *d_p* = (1/p)∏_{q<p}(1−1/q) and Σ_p *d_p* = 1). This code certifies the paper's three unconditional
perfect-packing results and generates its figures and empirical data:

* **Rectangle (§5).** The aspect-√2 rectangle is perfectly packed by similar rectangles of areas *d_p* —
  finite head certified in 200-bit interval arithmetic, infinite tail bounded by explicit prime estimates.
* **Ellipse (§6).** The ellipse of aspect ratio ρ ≈ 0.6296 is perfectly packed by similar ellipses of areas
  *d_p* in two orientations — a 1650-piece certified head (exact rational Perram–Wertheim disjointness
  witnesses, interval containment) plus the eroded-volume capacity tail, and the ten further certified
  *state heads* that carry the recursive integer dissection into elliptical cells.
* **Triangle (§7).** The isosceles right triangle is perfectly packed by similar triangles of areas *d_p*
  in the eight 45° orientations — a 10906-piece certified head (exact-rational and 100-bit interval
  containment margins, separating-axis disjointness witnesses on a rational-edge frame) plus the same
  eroded-volume capacity tail, now with the triangle's mixed-area constant κ_T = 3/√2, and the certified
  family of **45 per-state heads** (states 2…197) that carry the recursive integer dissection into
  triangular cells. An exact 15-piece algebraic base, verified symbolically in ℚ(√2,…,√47), anchors the
  head's first pieces.
* **Staircase (§3).** The staircase limit-angle constant L and θ\* = arctan(e^L) ≈ 56.29°.

---

## Contents

### Rectangle proof (§5)

| File | Role |
|------|------|
| `rect_exact.py` | The similar-rectangle model: guillotine scheduler, immature-stock measurement, and the two rigorous certificates behind the rectangle theorem (`head`, `tail`; the `tail` finite windows and analytic-tail constants run in `mpmath` interval arithmetic, and `head` derives its `CERTIFIED` status from the recorded invariants rather than asserting it). Also `ledger` and `backstop` float exploration runs. |
| `head_certificate_sqrt2.json` | The 5133-prime interval certificate for the √2 head, pre-generated (schema below). |

### Ellipse proof (§6)

| File | Role |
|------|------|
| `ellipse_certify.py` | The rigorous ellipse certifier. `certify` reads a head JSON and proves: area shares equal the exact σ_s(q) (reconstructed from *(state, n)*, never parsed from the file); every pair disjoint by an **exact rational** Perram–Wertheim witness F(λ\*) > 1 evaluated in `Fraction`; every piece contained by adaptive interval subdivision of the support gap. Fail-closed status. `tail` evaluates the placeability criterion, reports the crossover q₀ = 13967 and head size π(q₀) = 1650, and checks the analytic majorant. |
| `ellipse_build.py` | The head builder: beam search with spatial hashing, Sobol candidate fields, Nelder–Mead refinement gated on the exact margin the certifier will see, checkpoint/resume, and jam repair with escalating backup windows. `--state s` builds the state-*s* heads. Output JSON feeds `ellipse_certify.py` directly. |
| `ellipse_head_n1650.json` | The certified **container head** (state 2): 1650 pieces, primes 2…13967, exact dyadic centers, orientation pattern, ρ = 0.6296271949192772. |
| `head_s3.json` … `head_s31.json` | The ten further certified **state heads** of the state-comparison proposition (states 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 with n = 1044, 755, 594, 474, 405, 321, 273, 209, 171, 115). Every state s ≥ 37 is covered by the tail from its first piece, so these eleven heads are exactly what the recursive dissection of ℕ into elliptical cells requires. |

### Triangle proof (§7)

| File | Role |
|------|------|
| `triangle_certify.py` | The rigorous triangle certifier. `certify` reads a head JSON and proves, fail-closed: area shares equal the exact σ_s(q) (reconstructed from *(state, n)*, never parsed); every piece contained by **exact-rational** margins where a contact's radical coefficient vanishes (the intended leg/hypotenuse contacts) and 100-bit `mpmath` interval margins otherwise; every pair disjoint by a **separating-axis witness** among the eight directions *k*·45° (which contain every edge normal). `tail` evaluates the eroded-volume criterion with κ_T = 3/√2, reports the disciplined crossover q₀ = 115363 and head size π(q₀) = 10906, and checks the analytic majorant. |
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
  certification, and the triangle containment margins (`triangle_certify.py`, 100-bit).
* `sympy` — exact symbolic computation in ℚ(√2,…,√47) (`triangle_base_exact.py`) and exact prime/density
  bookkeeping for the ellipse and triangle scripts.
* `shapely` — floating-point polygon geometry (GEOS); used **only** to select packing combinatorics and to
  render figures. Every certified inequality is decided exactly (`sympy` / `Fraction`) or by
  directed-rounding interval arithmetic (`mpmath`), never by `shapely`.
* `scipy` — Sobol sampling and Nelder–Mead refinement in `ellipse_build.py` (search only; nothing certified
  depends on it).
* `numpy`, `matplotlib` — sieving, statistics, and figures.

Per script: `rect_exact.py` needs only `numpy` + `mpmath`; `ellipse_certify.py` needs `numpy` + `mpmath` +
`sympy`; `ellipse_build.py` needs `numpy` + `scipy` + `sympy`; `triangle_certify.py` needs `numpy` + `mpmath`
+ `sympy`; `triangle_build.py` needs only `numpy`; `run_triangle_states.py` needs `numpy` + `sympy` (it imports
`triangle_build` and shells out to `triangle_certify`); `triangle_base_exact.py` needs `numpy` + `sympy` +
`mpmath` + `shapely` (and its `_greedy.py` core); the figure pipeline needs `shapely` + `matplotlib` +
`sympy`.

---

## Quick start — reproduce the certificates

### Rectangle (§5)

```bash
# Finite head: place the first 5133 primes (q <= 5*10^4) in the sqrt2 rectangle,
# certify every placement fits in 200-bit interval arithmetic, and dump the full record.
python rect_exact.py head --q0 50000 --prec 200 --dump head_certificate_sqrt2.json

# Analytic tail: verify the immature-stock ratio rho(q) < 1 for all q > q0
# (finite windows in interval arithmetic + explicit analytic majorant; Appendix A).
python rect_exact.py tail
```

Expected head summary: `CERTIFIED`, `placements certified: 5133` (all primes 2…49999; piece 2 is the
implicit left-half split), reserve hit once (at *p* = 3) with slack 0.092, 5135 fringe rectangles, reserve
aspect √2 to the interval, reserve+fringe area = R_{<q₀} exactly, minimum fit margin 2.4×10⁻³.
Expected tail summary: `CERTIFIED`, ρ(q) ≤ 0.923 < 1−1/q for every prime q > q₀ = 5·10⁴.
Together these are the two halves of the unconditional, computer-assisted rectangle theorem.

### Ellipse (§6)

```bash
# Certify the 1650-piece container head: exact sigma_2(q) areas, exact rational
# Perram-Wertheim disjointness for every close pair, interval containment for every piece.
python ellipse_certify.py certify --in ellipse_head_n1650.json

# Certify any of the ten state heads (states 3..31) the same way:
python ellipse_certify.py certify --in head_s7.json

# Tail crossover: the placeability criterion holds for every prime q > q0 = 13967
# (head size pi(q0) = 1650), plus the analytic majorant checks.
python ellipse_certify.py tail
```

Expected certify output: `STATUS: CERTIFIED` with all pairs disjoint (exact PW witnesses for close pairs,
rigorous distance pruning for far pairs) and a certified negative max support gap. The container head
covers area Σ_{p≤13967} d_p = 1 − R_{<13997}; together with the `tail` crossover this is the
computer-assisted half of the ellipse theorem.

### Triangle (§7)

```bash
# Certify the 10906-piece container head: exact sigma_2(q) areas, exact-rational or
# 100-bit interval containment margins, separating-axis disjointness witnesses.
python triangle_certify.py certify --in triangle_head_n10906.json --design disciplined

# Tail crossover: the eroded-volume criterion (kappa_T = 3/sqrt2) holds for every
# prime q > q0 = 115363 (head size pi(q0) = 10906), plus the analytic majorant.
python triangle_certify.py tail --design disciplined

# Build + certify the whole 45-state integer family in one pass (writes states/).
# Drop the container head in first so state 2 is reused, not rebuilt:
cp triangle_head_n10906.json states/tri_state_2.json
python run_triangle_states.py --design disciplined --p0 1000 --outdir states
```

Expected certify output: `STATUS: CERTIFIED` with all 98154 containment margins ≥ 0 (11 exact-zero leg
contacts), all pairs disjoint by separating-axis witnesses, least certified gap ≈ 3.5×10⁻¹⁸. Expected family
output: `45/45 states CERTIFIED, 118593 pieces`, `FAMILY CERTIFIED`, least gap over the family ≈ 4.2×10⁻¹⁸.
Together with the `tail` crossover these are the computer-assisted half of the unconditional triangle theorem
and its integer corollary.

```bash
# Optional stricter check: the exact algebraic base (first 15 primes) in Q(sqrt2..sqrt47).
python triangle_base_exact.py                       # float greedy fixes the combinatorics, then exact proof
python triangle_base_exact.py --frozen              # replay from base_combinatorics.json (no greedy/GEOS)
```

Expected: `CERTIFIED` with `135/135`, `105/105`, free area = R_{<50}; `--frozen` gives the same result from
the stored record alone.

### Staircase constant (§3)

```bash
python staircase_constant.py                        # nmax = 10^8 (paper's bound); --nmax 10000000 is fast
```

Expected: `L = 0.4047892180`, `kappa = 1.4989865074`, `theta* = 56.2920568°`, `A = 0.2696063520`
(OEIS A078437), with the Leibniz bracket certifying `L ∈ [0.40478921, 0.40478923]`.

---

## `rect_exact.py` — the rectangle theorem

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

The certificate is a comparison of side **lengths**, so it is position-independent; a positioned rendering
of the same packing is the `rect_packing` figure.

---

## `ellipse_certify.py` and `ellipse_build.py` — the ellipse theorem

```
python ellipse_certify.py certify [--in ellipse_head_n1650.json] [--cells 96]
python ellipse_certify.py tail    [--rho 0.6296271949192772]
```

* **`certify`** — loads a head JSON, treating the stored binary64 centers as **exact dyadic rationals** and
  reconstructing every area share as the exact rational σ_s(q) = d_q/R_{<s} from *(state, n)* (their exact
  denominators are astronomical and are never parsed from the file). It then proves, fail-closed:
  * every close pair disjoint by an exact `Fraction` evaluation of the Perram–Wertheim contact function at
    a rational witness λ\* (far pairs pruned by a rigorous center-distance bound);
  * every piece contained by adaptive interval bisection (`mpmath.iv`) of the support gap in angle — a
    cell is accepted only when its interval upper bound is < 0.
* **`tail`** — the placeability criterion 2κ(1+T(q))/√(qR_{<q}) + N_q/q < 1 with κ = (ρ+1/ρ)/2; reports the
  crossover q₀ = 13967 (so the head is the first π(q₀) = 1650 primes) and evaluates the analytic majorant
  at 3·10⁴, 10⁵, 2·10⁵.

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
python triangle_certify.py tail    [--design disciplined] [--p0 1000] [--X 4000000] [--state 2]
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
* **`tail`** — the eroded-volume criterion (*b* + 2 Σ √d_p·V_p)/√(q·R_{<q}) + N_q/q < 1 with the mixed-area
  table *V_p* and κ_T = 3/√2; reports the disciplined crossover q₀ = 115363 (head π(q₀) = 10906) and checks
  the analytic majorant. `--state s` runs the criterion in a state-*s* cell (relative shares σ_s(q)).
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
| √2 rectangle packing (§5, Theorem `thm:rect`), finite head | `python rect_exact.py head --q0 50000 --dump head_certificate_sqrt2.json` |
| √2 rectangle packing, analytic tail (Lemma `lem:explicit` / App. A) | `python rect_exact.py tail` |
| Ellipse certified head (§6, Theorem `thm:ellhead`) | `python ellipse_certify.py certify --in ellipse_head_n1650.json` |
| Ellipse tail crossover (Lemma `lem:ellthresh`, Theorem `thm:ellipse`) | `python ellipse_certify.py tail` |
| Eleven state heads (§6, Prop. `prop:states`) | `python ellipse_certify.py certify --in head_sN.json` for N ∈ {3,5,7,11,13,17,19,23,29,31} |
| Integers in elliptical cells (§6, Theorem `thm:ellintegers`) | the eleven certified heads + `python ellipse_recursion_figure.py …` |
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

The certificates `head_certificate_sqrt2.json`, `ellipse_head_n1650.json` and the ten state heads
`head_s3.json`–`head_s31.json`, and `triangle_head_n10906.json` with the 45 state heads under `states/`
(and their `manifest.json`), together with the finite-window interval evaluations, are the archival record
behind the paper's three computer-assisted theorems.

## License

MIT — see `LICENSE`.
