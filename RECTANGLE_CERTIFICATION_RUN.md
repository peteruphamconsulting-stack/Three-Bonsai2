# Sharp LPF rectangle theorem: reproducible certification

These five scripts are required in the same directory:

1. `rectangle_full_interval_certificate.py` (main driver)
2. `rectangle_aspect_path_range.py` (interval-certified schedule ranges and algebraic endpoint contacts)
3. `rectangle_aspect_sweep.py` (prime enumeration and schedule support)
4. `rectangle_aspect_recursive_certificate.py` (shifted square-state verification)
5. `rectangle_aspect_head_certificate.py` (interval verification of state heads)

Prerequisite on a Mac with Python 3:

```bash
python3 -m pip install mpmath
```

Copy the five scripts into a common directory, open Terminal in that directory, and run:

```bash
python3 rectangle_full_interval_certificate.py
```

There are no required command-line flags. On success the final lines contain:

```
CLOSED-INTERVAL COVER: CERTIFIED (16 schedules plus two exact contacts)
UNIFORM MIXED-ORIENTATION TAIL: CERTIFIED
ALL SHIFTED STATES s>=3: CERTIFIED (universal aligned-square family)
OVERALL: CERTIFIED COMPLETE ASPECT INTERVAL
Saved rectangle_full_interval_certificate.json
```

The driver creates or replaces `rectangle_full_interval_certificate.json` in the same directory. The JSON records the rational inner intervals, schedule hashes, numerical tail bounds and status; it is accompanied by the code required to reconstruct and check the finite schedules. Keep all five Python files with the generated JSON in the research repository.

**Local cross-check:** The full driver was executed successfully on 2026-11-09. It verified all sixteen schedules, endpoint checks (2563 strict constraints each), the finite tail (maximum upper endpoint 0.9998320673924572 at prime 10487), the infinite majorant (0.7616246752362302), and the square-state family. This is verification by the supplied code, not an independent implementation of its mathematical logic.

Ellipse tail verification is separate. Use the revised `ellipse_certify.py` together with `revised_tail_certificate.py` in one directory:

```bash
python3 ellipse_certify.py tail
```

The ellipse command verifies the averaged mixed-area tail beyond the existing head through prime cutoff `10_000_000`, then an explicit infinite-range majorant. It does not regenerate geometric head files. Existing ellipse heads should still be checked using `python3 ellipse_certify.py certify --in ellipse_head_n1650.json` and analogous per-state files when available.
