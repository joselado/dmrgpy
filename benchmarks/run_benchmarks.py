#!/usr/bin/env python3
"""Cross-backend performance benchmark for DMRGPY.

Runs the same three calculations -- ground-state energy, a static
correlator, and a KPM dynamical correlator -- on a uniform S=1/2
Heisenberg chain, once per available backend (ITensor v2/v3 C++,
the pure-Python engine, and the live Julia backend when present),
over a sweep of chain lengths. Wall time is recorded for each
(backend, calculation, N) combination, together with an accuracy
cross-check against exact diagonalization at every size cheap enough
for ED, so that a large gap in *time* is never mistaken for a good
result if the backend is also silently computing the wrong physics.

The intent is to answer "which mode is worth optimizing next", not to
lock in golden regression values (see tests/test_benchmarks.py and the
various *_VS_ED examples for that) -- so a backend/size combination
that raises is recorded as a failure (ok=False) and the sweep
continues, rather than aborting the whole run.

Usage:
    python run_benchmarks.py
    python run_benchmarks.py --sizes 4 6 8 --dyn-max-n 6 --backends v3 python
    python run_benchmarks.py --repeats 3 --outdir /tmp/bench_out --no-pdf

Output: <outdir>/report.tex (and report.pdf if pdflatex is available),
<outdir>/results.json (raw data, so a report can be regenerated with
reporttk.write_report() without rerunning the sweep), and one plot per
calculation type under <outdir>/.
"""
import os
import sys
import time
import json
import shutil
import argparse
import datetime

# Insert this checkout's own src/ at the *front* of sys.path so this
# always exercises the current working tree, never a stale symlinked
# site-packages copy from a different checkout (see CLAUDE.md's
# "Worktree/symlink pitfall" -- same fix as tests/conftest.py).
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SRC = os.path.join(_ROOT, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
if os.path.dirname(os.path.abspath(__file__)) not in sys.path:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from dmrgpy import spinchain, cppext

import reporttk


# ---------------------------------------------------------------------
# Backend availability
# ---------------------------------------------------------------------

# (key, human label, itensor_version value accepted by Many_Body_Chain)
BACKEND_SPECS = [
    ("v2", "ITensor v2 (C++)", 2),
    ("v3", "ITensor v3 (C++)", 3),
    ("python", "Pure Python", "python"),
    ("julia_live", "Julia (ITensors.jl)", "julia_live"),
]

# julia_live is opt-in only: it's a fine backend once warmed up, but
# paying juliacall's per-process JIT compilation on every default run
# (even a run that never asked for Julia numbers) isn't worth it -- the
# caller must name it explicitly via --backends to get it included.
DEFAULT_BACKEND_KEYS = [k for k, _, _ in BACKEND_SPECS if k != "julia_live"]


def _julia_available():
    """Whether the live Julia backend can actually be used (the import
    itself is what boots the Julia session), mirroring cppext.available()
    for the compiled C++ backends -- same check used by
    examples/groundstate/dmrg_generalized_benchmark."""
    try:
        from dmrgpy.mpsjulialive import juliasession  # noqa: F401
        return True
    except Exception:
        return False


def available_backends(requested=None):
    """Return the BACKEND_SPECS rows that are actually usable in this
    environment, optionally restricted to a set of requested keys."""
    out = []
    for key, label, iv in BACKEND_SPECS:
        if requested is not None and key not in requested:
            continue
        if iv == 2 and not cppext.available(2):
            continue
        if iv == 3 and not cppext.available(3):
            continue
        if iv == "julia_live" and not _julia_available():
            continue
        out.append((key, label, iv))
    return out


# ---------------------------------------------------------------------
# Timing helper
# ---------------------------------------------------------------------

def timed(fn, repeats=1):
    """Run fn() up to `repeats` times, returning (value, best_time, ok,
    err). `value`/`ok` reflect the last successful call; `best_time` is
    the minimum wall time over all attempts (standard practice for wall
    time benchmarking -- robust to scheduler/GC noise, unlike a mean)."""
    value, err, ok, best = None, None, False, None
    for _ in range(max(1, repeats)):
        t0 = time.time()
        try:
            v = fn()
            dt = time.time() - t0
            value, ok, err = v, True, None
        except Exception as e:
            dt = time.time() - t0
            if not ok:
                err = repr(e)
        best = dt if best is None else min(best, dt)
    return value, best, ok, err


# ---------------------------------------------------------------------
# Problem construction
# ---------------------------------------------------------------------

def make_chain(n, itensor_version, maxm, nsweeps, cutoff):
    """A uniform S=1/2 Heisenberg chain of length n on the requested
    backend -- the standard test Hamiltonian used across this codebase's
    own examples/tests (see e.g. test_benchmarks.py's Heisenberg test)."""
    spins = ["1/2" for _ in range(n)]
    sc = spinchain.Spin_Chain(spins, itensor_version=itensor_version)

    def fj(i, j):
        return 1.0 if abs(i - j) == 1 else 0.0

    sc.set_exchange(fj)
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, nsweeps, cutoff
    return sc


def ed_reference(n, params):
    """Ground-state energy and <Sz_0 Sz_i> static correlator from exact
    diagonalization -- used only as an accuracy cross-check, independent
    of itensor_version (any chain object's ED path is identical), so the
    pure-Python backend is used here since it needs no compiled
    extension to even construct the chain."""
    sc = make_chain(n, "python", params["maxm"], params["nsweeps"], params["cutoff"])
    e, _, ok_e, err_e = timed(lambda: sc.gs_energy(mode="ED"))
    if not ok_e:
        return dict(n=n, ok=False, err=err_e, energy=None, correlator=None)
    corr, _, ok_c, err_c = timed(
        lambda: np.array([sc.vev(sc.Sz[0] * sc.Sz[i], mode="ED").real for i in range(n)]))
    return dict(n=n, ok=ok_c, err=err_c, energy=e,
                correlator=corr if ok_c else None)


# ---------------------------------------------------------------------
# One (backend, N) case: ground state + static correlator (+ dynamical
# correlator when requested), all sharing the one DMRG ground state.
# ---------------------------------------------------------------------

def run_case(key, iv, n, params, run_dynamic):
    rows = []
    repeats = params["repeats"]

    sc, _, ok_build, err_build = timed(
        lambda: make_chain(n, iv, params["maxm"], params["nsweeps"], params["cutoff"]))
    if not ok_build:
        rows.append(dict(calc="gs", backend=key, n=n, ok=False, err=err_build,
                          time=None, energy=None))
        rows.append(dict(calc="static", backend=key, n=n, ok=False,
                          err="chain construction failed", time=None, value=None))
        if run_dynamic:
            rows.append(dict(calc="dynamic", backend=key, n=n, ok=False,
                              err="chain construction failed", time=None, value=None))
        return rows

    e, t_gs, ok_gs, err_gs = timed(lambda: sc.gs_energy(mode="DMRG"), repeats=repeats)
    rows.append(dict(calc="gs", backend=key, n=n, ok=ok_gs, err=err_gs, time=t_gs, energy=e))
    if not ok_gs:
        rows.append(dict(calc="static", backend=key, n=n, ok=False,
                          err="ground state failed", time=None, value=None))
        if run_dynamic:
            rows.append(dict(calc="dynamic", backend=key, n=n, ok=False,
                              err="ground state failed", time=None, value=None))
        return rows

    # Static correlator: <Sz_0 Sz_i> for every site, computed on top of
    # the already-converged (and now cached) ground state above, so this
    # timing reflects the correlator cost alone, not gs_energy's cost.
    corr, t_static, ok_s, err_s = timed(
        lambda: np.array([sc.vev(sc.Sz[0] * sc.Sz[i]).real for i in range(n)]),
        repeats=repeats)
    rows.append(dict(calc="static", backend=key, n=n, ok=ok_s, err=err_s,
                      time=t_static, value=corr))

    if run_dynamic:
        dc, t_dyn, ok_d, err_d = timed(
            lambda: sc.get_dynamical_correlator(
                mode="DMRG", submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                n=params["kpm_n"], es=params["kpm_es"]),
            repeats=repeats)
        rows.append(dict(calc="dynamic", backend=key, n=n, ok=ok_d, err=err_d,
                          time=t_dyn, value=(dc[1].real if ok_d else None)))
    return rows


# ---------------------------------------------------------------------
# CLI / orchestration
# ---------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[6, 10, 14, 18, 22, 26, 30],
            help="chain lengths swept for gs_energy and the static correlator")
    ap.add_argument("--dyn-max-n", type=int, default=20,
            help="largest N (from --sizes) also run through the KPM "
                 "dynamical correlator, which is the most expensive of "
                 "the three calculations")
    ap.add_argument("--ed-max-n", type=int, default=18,
            help="largest N for which the exact-diagonalization accuracy "
                 "cross-check is even attempted -- ED's 2**N-dimensional "
                 "Hilbert space makes it infeasible well before DMRG/KPM "
                 "stop being cheap, so sizes above this are timed but not "
                 "accuracy-checked")
    ap.add_argument("--backends", nargs="+", default=None,
            choices=[k for k, _, _ in BACKEND_SPECS],
            help="restrict to these backend keys (default: every backend "
                 "available in this environment EXCEPT julia_live, which "
                 "is opt-in only -- pass it explicitly, e.g. "
                 "'--backends v2 v3 python julia_live', to include it; "
                 "juliacall's per-process JIT warm-up cost makes it not "
                 "worth paying by default on every run)")
    ap.add_argument("--maxm", type=int, default=40)
    ap.add_argument("--nsweeps", type=int, default=8)
    ap.add_argument("--cutoff", type=float, default=1e-8)
    ap.add_argument("--kpm-n", type=int, default=300, help="KPM moment count")
    ap.add_argument("--kpm-emin", type=float, default=-1.0)
    ap.add_argument("--kpm-emax", type=float, default=6.0)
    ap.add_argument("--kpm-npts", type=int, default=120)
    ap.add_argument("--repeats", type=int, default=1,
            help="repeat each timed call and keep the minimum wall time")
    ap.add_argument("--outdir", default=None,
            help="output directory (default: benchmarks/output next to this script)")
    ap.add_argument("--no-pdf", action="store_true",
            help="skip compiling report.tex with pdflatex even if it's on PATH")
    ap.add_argument("--from-json", default=None,
            help="skip the sweep entirely and regenerate report.tex from a "
                 "previously written results.json (e.g. after tweaking "
                 "reporttk.py) -- all other sweep options are ignored")
    return ap.parse_args()


def _rerender_from_json(json_path, outdir, no_pdf):
    with open(json_path) as f:
        data = json.load(f)
    backends = [(k, l, iv) for (k, l) in data["backends"]
            for _k2, _l2, iv in BACKEND_SPECS if _k2 == k]
    params = dict(data["params"])
    params["kpm_es"] = np.array(params["kpm_es"])
    ed_refs = {int(n): ref for n, ref in data["ed_refs"].items()}
    tex_path = reporttk.write_report(data["results"], ed_refs, params, backends,
            data["sizes"], data["dyn_sizes"], outdir,
            generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print(f"Wrote {tex_path}")
    if not no_pdf and shutil.which("pdflatex"):
        ok, msg = reporttk.compile_pdf(tex_path)
        print(msg if ok else f"pdflatex failed:\n{msg}")


def main():
    args = parse_args()
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(outdir, exist_ok=True)

    if args.from_json:
        _rerender_from_json(args.from_json, outdir, args.no_pdf)
        return

    params = dict(maxm=args.maxm, nsweeps=args.nsweeps, cutoff=args.cutoff,
            kpm_n=args.kpm_n,
            kpm_es=np.linspace(args.kpm_emin, args.kpm_emax, args.kpm_npts),
            repeats=args.repeats)

    requested = args.backends if args.backends is not None else DEFAULT_BACKEND_KEYS
    backends = available_backends(requested=requested)
    if not backends:
        print("No requested backend is available in this environment -- nothing to do.")
        sys.exit(1)

    requested_but_missing = set(requested) - {k for k, _, _ in backends}
    for key in sorted(requested_but_missing):
        print(f"(backend {key!r} was requested but is not available here -- skipping)")
    if args.backends is None:
        print("(julia_live skipped by default -- pass '--backends ... julia_live' to include it)")

    sizes = sorted(set(args.sizes))
    dyn_sizes = set(n for n in sizes if n <= args.dyn_max_n)

    print("Backends:  " + ", ".join(f"{label} ({key})" for key, label, _ in backends))
    print("Sizes:     " + ", ".join(str(n) for n in sizes))
    print("Dynamical: " + ", ".join(str(n) for n in sorted(dyn_sizes)) +
          "  (KPM correlator, the most expensive calculation)")
    print()

    if any(key == "julia_live" for key, _, _ in backends):
        # juliacall JIT-compiles each distinct Julia method signature the
        # first time it's called in a process -- confirmed in this
        # codebase to cost tens of seconds, independent of the problem
        # size used to trigger it (see CLAUDE.md's "julia_live tests only
        # need to run..." section). Without this warm-up, the very first
        # julia_live timing in the sweep below would be dominated by JIT
        # compilation rather than the algorithm, and would badly skew the
        # relative-speed summary. One throwaway tiny run before the timed
        # sweep starts pays that one-time cost up front instead.
        print("Warming up the Julia backend (one-time JIT compilation, not timed)...")
        try:
            run_case("julia_live", "julia_live", min(sizes), params, bool(dyn_sizes))
        except Exception as e:
            print(f"  (warm-up run raised {e!r} -- continuing anyway)")
        print()

    results = []
    t_start = time.time()
    for key, label, iv in backends:
        for n in sizes:
            run_dynamic = n in dyn_sizes
            print(f"[{time.time()-t_start:7.1f}s] {label:22s} N={n:3d} ...", end="", flush=True)
            rows = run_case(key, iv, n, params, run_dynamic)
            results.extend(rows)
            by_calc = {r["calc"]: r for r in rows}
            parts = []
            for calc in ("gs", "static", "dynamic"):
                if calc not in by_calc:
                    continue
                r = by_calc[calc]
                parts.append(f"{calc}={r['time']:.3f}s" if r["ok"] else f"{calc}=FAIL")
            print("  " + "  ".join(parts))

    print(f"\nTotal sweep time: {time.time()-t_start:.1f}s")

    print("\nComputing ED reference values for the accuracy cross-check "
          f"(N <= {args.ed_max_n})...")
    ed_refs = {}
    for n in sizes:
        if n > args.ed_max_n:
            print(f"  N={n:3d}: skipped (above --ed-max-n={args.ed_max_n}, "
                  f"2**{n} is infeasible for exact diagonalization)")
            continue
        ref = ed_reference(n, params)
        ed_refs[n] = ref
        status = "ok" if ref["ok"] else f"FAILED ({ref['err']})"
        print(f"  N={n:3d}: {status}")

    with open(os.path.join(outdir, "results.json"), "w") as f:
        json.dump(reporttk.to_jsonable(dict(
            results=results, ed_refs=ed_refs, params=params,
            backends=[(k, l) for k, l, _ in backends], sizes=sizes,
            dyn_sizes=sorted(dyn_sizes))), f, indent=2)
    print(f"Wrote {os.path.join(outdir, 'results.json')}")

    tex_path = reporttk.write_report(results, ed_refs, params, backends,
            sizes, sorted(dyn_sizes), outdir,
            generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print(f"Wrote {tex_path}")

    if not args.no_pdf and shutil.which("pdflatex"):
        ok, msg = reporttk.compile_pdf(tex_path)
        print(msg if ok else f"pdflatex failed:\n{msg}")
    elif not args.no_pdf:
        print("(pdflatex not found on PATH -- leaving report.tex uncompiled)")


if __name__ == "__main__":
    main()
