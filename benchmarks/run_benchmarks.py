#!/usr/bin/env python3
"""Cross-backend performance benchmark for DMRGPY.

Runs the same set of calculations on the uniform S=1/2 Heisenberg chain,
once per available backend (ITensor v2/v3 C++, the pure-Python engine,
and the live Julia backend when present), and records wall time for every
(backend, calculation, size) combination together with an accuracy
cross-check, so that a large gap in *time* is never mistaken for a good
result if the backend is also silently computing the wrong physics.

Two families of calculation are swept (see calctk.py, which is the single
list of what is covered):

* **finite chains**, swept over chain length N -- ground-state energy,
  static correlator, KPM and TD dynamical correlators, first excited
  state, bond entanglement entropy, real-time evolution after a quench,
  and the ground state confined to the Sz=0 conserved sector. Accuracy is
  checked against exact diagonalization at every size cheap enough for it.
* **infinite chains**, swept over bond dimension chi -- energy density
  and <Sz> under both gs_method="idmrg" and gs_method="vumps". Accuracy is
  checked against the Bethe-ansatz energy density and against <Sz>=0.

The intent is to answer "which mode is worth optimizing next", not to
lock in golden regression values (see tests/test_benchmarks.py and the
various *_VS_ED examples for that) -- so a backend/size combination
that raises is recorded as a failure (ok=False) and the sweep
continues, rather than aborting the whole run. Calculations that only
some backends implement (the conserved sector, the infinite chain) are
part of the sweep too, restricted to the backends that have them.

Usage:
    python run_benchmarks.py
    python run_benchmarks.py --sizes 4 6 8 --dyn-max-n 6 --backends v3 python
    python run_benchmarks.py --calcs gs static --repeats 3 --no-pdf

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

from dmrgpy import cppext

import calctk
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

def timed(make_thunk, repeats=1):
    """Time a calculation up to `repeats` times, returning (value,
    best_time, ok, err). `value`/`ok` reflect the last successful call;
    `best_time` is the minimum wall time over all attempts (standard
    practice for wall-time benchmarking -- robust to scheduler/GC noise,
    unlike a mean).

    `make_thunk(rep)` is called once per repeat and does whatever untimed
    preparation that repeat needs, returning the callable that is
    actually timed. It gets the repeat index because some calculations
    are cached: the DMRG session caches its ground-state energy
    (groundstate.py's send-cache), so a second gs_energy() on the same
    chain returns instantly, and timing that would report a backend as
    infinitely fast. Preparation per repeat lets those calculations hand
    back a fresh chain instead."""
    value, err, ok, best = None, None, False, None
    for rep in range(max(1, repeats)):
        try:
            fn = make_thunk(rep)
        except Exception as e:
            if not ok:
                err = repr(e)
            continue
        t0 = time.time()
        try:
            v = fn()
            dt = time.time() - t0
            value, ok, err = v, True, None
        except Exception as e:
            dt = time.time() - t0
            if not ok:
                err = repr(e)
            continue
        best = dt if best is None else min(best, dt)
    return value, best, ok, err


def _fail_row(calc, key, n, err):
    return dict(calc=calc, backend=key, n=n, ok=False, err=err, time=None,
                value=None)


# ---------------------------------------------------------------------
# One (backend, N) case of the finite family: every selected finite
# calculation, sharing the one DMRG ground state where they can.
# ---------------------------------------------------------------------

# Calculations measured on top of the shared, already-converged ground
# state of this (backend, N) -- as opposed to the ones that build their
# own chain because they change the Hamiltonian (evolution) or the site
# set (sector).
_NEEDS_GS = ("static", "dynamic", "dynamic_td", "excited", "entropy")


def run_case_finite(key, iv, n, params, calcs):
    rows = []
    if not calcs:
        return rows

    try:
        sc = calctk.make_chain(n, iv, params)
    except Exception as e:
        return [_fail_row(c, key, n, "chain construction failed: %r" % (e,))
                for c in calcs]

    gs_ok = None  # None = not attempted yet
    for calc in calcs:
        prep = calctk.spec(calc)["prep"]
        if calc in _NEEDS_GS:
            if gs_ok is None:  # "gs" itself wasn't part of this run
                _e, _t, gs_ok, _err = timed(
                        lambda rep: (lambda: sc.gs_energy(mode="DMRG")))
            if not gs_ok:
                rows.append(_fail_row(calc, key, n, "ground state failed"))
                continue
        value, t, ok, err = timed(
                lambda rep, _p=prep: _p(sc, iv, n, params, rep),
                repeats=params["repeats"])
        if calc == "gs":
            gs_ok = ok
        rows.append(dict(calc=calc, backend=key, n=n, ok=ok, err=err,
                         time=t, value=value))
    return rows


# ---------------------------------------------------------------------
# One (backend, chi) case of the infinite family. The two gs_methods get
# one converged chain each, shared by that method's own observables.
# ---------------------------------------------------------------------

# The infinite family's own ground-state entries: everything else in
# that family is an observable measured on top of one of these.
_INFINITE_GS = ("idmrg_gs", "vumps_gs")


def run_case_infinite(key, iv, chi, params, calcs):
    rows = []
    for method in ("idmrg", "vumps"):
        todo = [c for c in calcs if calctk.spec(c).get("method") == method]
        if not todo:
            continue
        try:
            ic = calctk.make_infinite_chain(chi, iv, method, params)
        except Exception as e:
            rows += [_fail_row(c, key, chi, "chain construction failed: %r" % (e,))
                     for c in todo]
            continue
        gs_ok = None
        for calc in todo:
            prep = calctk.spec(calc)["prep"]
            if calc not in _INFINITE_GS:
                if gs_ok is None:  # the ground state wasn't itself timed
                    _e, _t, gs_ok, _err = timed(
                            lambda rep: (lambda: ic.gs_energy()))
                if not gs_ok:
                    rows.append(_fail_row(calc, key, chi, "ground state failed"))
                    continue
            value, t, ok, err = timed(
                    lambda rep, _p=prep: _p(ic, chi, iv, params, rep),
                    repeats=params["repeats"])
            if calc in _INFINITE_GS:
                gs_ok = ok
            rows.append(dict(calc=calc, backend=key, n=chi, ok=ok, err=err,
                             time=t, value=value))
    return rows


# ---------------------------------------------------------------------
# CLI / orchestration
# ---------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[6, 10, 14, 18, 22, 26, 30],
            help="chain lengths swept by the finite-chain calculations")
    ap.add_argument("--chis", type=int, nargs="+", default=[8, 16, 24, 32],
            help="bond dimensions swept by the infinite-chain calculations "
                 "(their sweep axis: an infinite chain has no length)")
    ap.add_argument("--calcs", nargs="+", default=None,
            choices=calctk.CALC_KEYS,
            help="restrict to these calculations (default: all of them; "
                 "see calctk.py for what each one measures)")
    ap.add_argument("--dyn-max-n", type=int, default=20,
            help="largest N run through the KPM dynamical correlator")
    ap.add_argument("--td-max-n", type=int, default=14,
            help="largest N run through the TD (real-time + Fourier) "
                 "dynamical correlator, the most expensive calculation here")
    ap.add_argument("--evolution-max-n", type=int, default=20,
            help="largest N run through the real-time quench evolution")
    ap.add_argument("--excited-max-n", type=int, default=22,
            help="largest N run through the excited-state calculation")
    ap.add_argument("--vumps-max-chi", type=int, default=24,
            help="largest bond dimension run through the VUMPS solver, "
                 "which grows out of proportion to everything else here")
    ap.add_argument("--ed-max-n", type=int, default=14,
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
    ap.add_argument("--emin", type=float, default=-1.0,
            help="lower edge of the energy window of the spectral functions")
    ap.add_argument("--emax", type=float, default=6.0,
            help="upper edge of that window")
    ap.add_argument("--npts", type=int, default=120,
            help="number of energies sampled in that window")
    ap.add_argument("--td-delta", type=float, default=0.2,
            help="Lorentzian broadening of the TD dynamical correlator, "
                 "which also sets how long it has to evolve for")
    ap.add_argument("--nt", type=int, default=100,
            help="time steps of the real-time evolution benchmark")
    ap.add_argument("--dt", type=float, default=0.05,
            help="time step of the real-time evolution benchmark")
    ap.add_argument("--inf-maxiter", type=int, default=60,
            help="iteration cap of the infinite-chain ground-state solvers")
    ap.add_argument("--inf-etol", type=float, default=1e-10,
            help="energy tolerance of the infinite-chain ground-state solvers")
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
    missing = [k for k in ("results", "refs", "calc_sizes", "ref_kinds")
               if k not in data]
    if missing:
        sys.exit("%s is missing %s -- it predates the current results.json "
                 "layout, so rerun the sweep instead of re-rendering it"
                 % (json_path, ", ".join(missing)))
    backends = [(k, l, iv) for (k, l) in data["backends"]
            for _k2, _l2, iv in BACKEND_SPECS if _k2 == k]
    params = dict(data["params"])
    params["es"] = np.array(params["es"])
    refs = {calc: {int(x): v for x, v in per_x.items()}
            for calc, per_x in data["refs"].items()}
    tex_path = reporttk.write_report(data["results"], refs, params, backends,
            data["calc_sizes"], data["ref_kinds"], outdir,
            generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print(f"Wrote {tex_path}")
    if not no_pdf and shutil.which("pdflatex"):
        ok, msg = reporttk.compile_pdf(tex_path)
        print(msg if ok else f"pdflatex failed:\n{msg}")


def _calc_sizes(args, calcs, sizes, chis):
    """The sweep axis of every selected calculation: chain lengths for
    the finite family, bond dimensions for the infinite one, each capped
    per calculation since they are nowhere near equally expensive (the
    real-time correlator and VUMPS in particular would otherwise
    dominate the wall time of a default run)."""
    caps = {"dynamic": args.dyn_max_n, "dynamic_td": args.td_max_n,
            "evolution": args.evolution_max_n, "excited": args.excited_max_n,
            "vumps_gs": args.vumps_max_chi, "vumps_vev": args.vumps_max_chi}
    out = {}
    for calc in calcs:
        axis = chis if calctk.spec(calc)["group"] == "infinite" else sizes
        cap = caps.get(calc)
        out[calc] = [x for x in axis if cap is None or x <= cap]
    return out


def _references(calc_sizes, params, calcs, ed_max_n):
    """{calc: {size: reference value}} for every calculation that has
    one. The finite family shares a single ED run per chain length; the
    infinite family's references are closed-form and size-independent."""
    refs = {c: {} for c in calcs if calctk.spec(c)["ref"] is not None}
    finite = [c for c in calcs if calctk.spec(c)["group"] == "finite"
              and calctk.spec(c)["ref"] is not None]
    ed_sizes = sorted(set(n for c in finite for n in calc_sizes[c] if n <= ed_max_n))
    for n in ed_sizes:
        ref = calctk.finite_references(n, params, finite)
        status = "ok" if ref["ok"] else "FAILED (%s)" % ref["err"]
        print(f"  N={n:3d}: {status}")
        if not ref["ok"]:
            continue
        for calc in finite:
            key = calctk.spec(calc)["ref"]
            if key in ref and n in calc_sizes[calc]:
                refs[calc][n] = ref[key]
    for calc in calcs:
        s = calctk.spec(calc)
        if s["group"] != "infinite" or s["ref"] is None:
            continue
        value = calctk.infinite_reference(s["ref"])
        for chi in calc_sizes[calc]:
            refs[calc][chi] = value
    skipped = [n for c in finite for n in calc_sizes[c] if n > ed_max_n]
    if skipped:
        print(f"  N > {ed_max_n} skipped (2**N is infeasible for exact "
              f"diagonalization): " + ", ".join(str(n) for n in sorted(set(skipped))))
    return refs


def main():
    args = parse_args()
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(outdir, exist_ok=True)

    if args.from_json:
        _rerender_from_json(args.from_json, outdir, args.no_pdf)
        return

    params = dict(maxm=args.maxm, nsweeps=args.nsweeps, cutoff=args.cutoff,
            kpm_n=args.kpm_n,
            es=np.linspace(args.emin, args.emax, args.npts),
            td_delta=args.td_delta, nt=args.nt, dt=args.dt,
            inf_maxiter=args.inf_maxiter, inf_etol=args.inf_etol,
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

    calcs = args.calcs if args.calcs is not None else calctk.CALC_KEYS
    calcs = [c for c in calctk.CALC_KEYS if c in calcs]  # canonical order
    sizes = sorted(set(args.sizes))
    chis = sorted(set(args.chis))
    calc_sizes = _calc_sizes(args, calcs, sizes, chis)

    print("Backends:  " + ", ".join(f"{label} ({key})" for key, label, _ in backends))
    for calc in calcs:
        s = calctk.spec(calc)
        axis = "chi" if s["group"] == "infinite" else "N"
        only = "" if s["backends"] is None else \
                "   [%s only]" % ", ".join(s["backends"])
        print(f"  {calc:12s} {axis}=" +
              ",".join(str(x) for x in calc_sizes[calc]) + only)
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
        warm = [c for c in calcs if calctk.spec(c)["group"] == "finite"
                and calctk.allows(c, "julia_live")]
        try:
            run_case_finite("julia_live", "julia_live", min(sizes), params, warm)
        except Exception as e:
            print(f"  (warm-up run raised {e!r} -- continuing anyway)")
        print()

    results = []
    t_start = time.time()
    for key, label, iv in backends:
        finite = [c for c in calcs if calctk.spec(c)["group"] == "finite"
                  and calctk.allows(c, key)]
        infinite = [c for c in calcs if calctk.spec(c)["group"] == "infinite"
                    and calctk.allows(c, key)]
        for n in sizes:
            todo = [c for c in finite if n in calc_sizes[c]]
            if not todo:
                continue
            print(f"[{time.time()-t_start:7.1f}s] {label:22s} N={n:3d} ...",
                  end="", flush=True)
            rows = run_case_finite(key, iv, n, params, todo)
            results.extend(rows)
            print("  " + "  ".join(
                (f"{r['calc']}={r['time']:.3f}s" if r["ok"] else f"{r['calc']}=FAIL")
                for r in rows))
        for chi in chis:
            todo = [c for c in infinite if chi in calc_sizes[c]]
            if not todo:
                continue
            print(f"[{time.time()-t_start:7.1f}s] {label:22s} chi={chi:3d} ...",
                  end="", flush=True)
            rows = run_case_infinite(key, iv, chi, params, todo)
            results.extend(rows)
            print("  " + "  ".join(
                (f"{r['calc']}={r['time']:.3f}s" if r["ok"] else f"{r['calc']}=FAIL")
                for r in rows))

    print(f"\nTotal sweep time: {time.time()-t_start:.1f}s")

    print("\nComputing reference values for the accuracy cross-check "
          f"(exact diagonalization for N <= {args.ed_max_n})...")
    refs = _references(calc_sizes, params, calcs, args.ed_max_n)
    # What each calculation's accuracy column is measured against, for
    # the report (reporttk.REF_HEADERS names them).
    ref_kinds = {c: calctk.spec(c)["ref"] for c in calcs
                 if calctk.spec(c)["ref"] is not None}

    with open(os.path.join(outdir, "results.json"), "w") as f:
        json.dump(reporttk.to_jsonable(dict(
            results=results, refs=refs, params=params,
            backends=[(k, l) for k, l, _ in backends],
            calc_sizes=calc_sizes, ref_kinds=ref_kinds)), f, indent=2)
    print(f"Wrote {os.path.join(outdir, 'results.json')}")

    tex_path = reporttk.write_report(results, refs, params, backends,
            calc_sizes, ref_kinds, outdir,
            generated_at=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print(f"Wrote {tex_path}")

    if not args.no_pdf and shutil.which("pdflatex"):
        ok, msg = reporttk.compile_pdf(tex_path)
        print(msg if ok else f"pdflatex failed:\n{msg}")
    elif not args.no_pdf:
        print("(pdflatex not found on PATH -- leaving report.tex uncompiled)")


if __name__ == "__main__":
    main()
