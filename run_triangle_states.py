#!/usr/bin/env python3
"""
run_triangle_states.py -- build and certify the full family of triangle state heads
(the sub-dissection heads that turn the triangle packing into a geometric SIEVE,
exactly as the eleven ellipse state heads do for the ellipse).

For each prime state s from 2 up to the terminal state (the first prime whose head
is empty), this builds the state-s head with triangle_build.py and certifies it with
triangle_certify.py, writing one head JSON per state plus a manifest. For the
disciplined design the states are the primes 2..197 (45 heads, ~118.6k pieces; the
base state 2 with 10906 pieces is the largest and every later head is smaller).

An existing CERTIFIED head at the target path is detected and reused, so an
interrupted run resumes cheaply and a pre-built base head (state 2) is not rebuilt:
just drop it at <outdir>/tri_state_2.json first.

Usage:
    python3 run_triangle_states.py --design disciplined --p0 1000 --outdir states
    python3 run_triangle_states.py --only 197,193,191       # a subset, for testing
    python3 run_triangle_states.py --keep-going             # don't stop on a failure
    python3 run_triangle_states.py --certify-only           # re-verify existing heads

Requires: numpy, sympy, mpmath (via the two scripts).
"""
import argparse, json, os, re, subprocess, sys, time
import numpy as np, sympy
import triangle_build as tb

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "triangle_build.py")
CERT  = os.path.join(HERE, "triangle_certify.py")

def run(cmd):
    return subprocess.run([sys.executable] + cmd, capture_output=True, text=True)

def certify(headfile, prec=100):
    r = run([CERT, "certify", "--in", headfile, "--prec", str(prec)])
    out = r.stdout + r.stderr
    status = "CERTIFIED" if "STATUS: CERTIFIED" in out else "FAILED"
    m = {}
    for key, pat in [("cont_min", r"min margin ([+\-0-9.eE]+)"),
                     ("sat_min",  r"min witness gap ([0-9.eE\-]+)"),
                     ("n",        r"n=(\d+)"), ("qlast", r"\.\.(\d+),")]:
        g = re.search(pat, out)
        if g: m[key] = g.group(1)
    return status, m, out

def build(state, design, p0, headfile, seed):
    ck = headfile + ".ckpt"
    cmd = [BUILD, "--state", str(state), "--design", design, "--p0", str(p0),
           "--out", headfile, "--checkpoint", ck, "--ckpt-every", "1000", "--seed", str(seed)]
    if os.path.exists(ck):
        cmd += ["--resume", ck]
    r = run(cmd); return r.stdout + r.stderr

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", default="disciplined", choices=["uniform8", "disciplined"])
    ap.add_argument("--p0", type=int, default=1000)
    ap.add_argument("--outdir", default="states")
    ap.add_argument("--only", default="")            # comma list of states
    ap.add_argument("--max-state", type=int, default=400)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--keep-going", action="store_true")
    ap.add_argument("--certify-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    P = tb.sieve(800_000)
    if a.only:
        states = [int(x) for x in a.only.split(",") if x.strip()]
    else:
        states = [int(s) for s in sympy.primerange(2, a.max_state+1)
                  if tb.head_count(P, int(s), a.design, a.p0) > 0]
    total = sum(tb.head_count(P, s, a.design, a.p0) for s in states)
    print(f"design={a.design} p0={a.p0}: {len(states)} state heads, "
          f"{total} pieces total, states {states[0]}..{states[-1]}")

    manifest = {"design": a.design, "p0": a.p0, "states": [], "total_pieces": total,
                "terminal_state": int(sympy.nextprime(states[-1]))}
    t0 = time.time(); all_ok = True
    for s in states:
        hf = os.path.join(a.outdir, f"tri_state_{s}.json")
        want = tb.head_count(P, s, a.design, a.p0)
        # reuse an already-certified head
        reused = False
        if os.path.exists(hf):
            st, m, _ = certify(hf)
            if st == "CERTIFIED" and int(m.get("n", -1)) == want:
                reused = True
                print(f"  s={s:>4}: reuse CERTIFIED head (n={want})")
        if not reused:
            if a.certify_only:
                print(f"  s={s:>4}: MISSING head, skipped (certify-only)"); all_ok = False; continue
            tb_t = time.time()
            build(s, a.design, a.p0, hf, a.seed)
            st, m, out = certify(hf)
            dt = time.time() - tb_t
            print(f"  s={s:>4}: built n={m.get('n','?')} in {dt:.0f}s -> {st}"
                  f"  (min margin {m.get('cont_min','?')}, min gap {m.get('sat_min','?')})")
            if st != "CERTIFIED":
                all_ok = False
                print(f"        !! {st}; see head {hf}")
                if not a.keep_going:
                    print("stopping (use --keep-going to continue)"); break
        manifest["states"].append({"state": s, "n": want, "head": os.path.basename(hf),
                                    "status": "CERTIFIED" if reused else st})
    manifest["all_certified"] = all_ok
    manifest["elapsed_sec"] = round(time.time() - t0, 1)
    mf = os.path.join(a.outdir, "manifest.json")
    json.dump(manifest, open(mf, "w"), indent=1)
    ncert = sum(1 for e in manifest["states"] if e["status"] == "CERTIFIED")
    print(f"\n{ncert}/{len(states)} states CERTIFIED, {time.time()-t0:.0f}s total; "
          f"manifest -> {mf}")
    print("FAMILY CERTIFIED" if all_ok and ncert == len(states) else "INCOMPLETE / FAILURES")

if __name__ == "__main__":
    main()
