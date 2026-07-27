# Head-to-head benchmark: arnolditk's thick-restart Arnoldi (the
# pre-existing MPS Arnoldi route, see dmrgpy.algebra.arnolditk) vs the new
# matrix-free Implicitly Restarted Arnoldi Method (IRAM) ported from
# ARPACK's znaupd/znaup2 (see dmrgpy.algebra.arpacktk), focusing on
# non-Hermitian ground states first (mode="GS" / which="SR" -- smallest
# real part) since that is the case arnolditk's own docs call out as
# "typically several orders of magnitude less accurate than NH-DMRG at
# comparable cost" (see examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi).
#
# Instruments the number of expensive Op(x) calls (MPO*MPS application in
# DMRG mode / sparse matvec in ED mode) alongside wall time and accuracy
# against an ED ground truth -- the same instrumentation approach as
# examples/groundstate/arnoldi_benchmark, extended to cover both
# algorithms so they can be compared directly on identical test cases.
#
# Usage: python main.py [--json out.json] [--modes ED,DMRG]
import sys, os, time, json, argparse

sys.path.append(os.getcwd()+'/../../../src')

import numpy as np

from dmrgpy.multioperatortk.staticoperator import StaticOperator
from dmrgpy.edtk.edchain import EDOperator, State
from dmrgpy.mps import MPS
from dmrgpy.algebra import arnolditk
from dmrgpy.algebra import arpacktk

# ---------------------------------------------------------------------
# Instrumentation: count Op(x) = M*x calls, the expensive primitive.
# ---------------------------------------------------------------------
_counts = {"dmrg": 0, "ed": 0}

_orig_static_mul = StaticOperator.__mul__
def _counted_static_mul(self, v):
    if type(v) == MPS:
        _counts["dmrg"] += 1
    return _orig_static_mul(self, v)
StaticOperator.__mul__ = _counted_static_mul

_orig_ed_mul = EDOperator.__mul__
def _counted_ed_mul(self, v):
    if isinstance(v, State):
        _counts["ed"] += 1
    return _orig_ed_mul(self, v)
EDOperator.__mul__ = _counted_ed_mul

def reset_counts():
    _counts["dmrg"] = 0
    _counts["ed"] = 0

def get_count(mode):
    return _counts["dmrg"] if mode == "DMRG" else _counts["ed"]


# ---------------------------------------------------------------------
# Test cases: non-Hermitian first (the focus of this comparison), then a
# couple of Hermitian ones since IRAM/arnolditk should agree there too.
# ---------------------------------------------------------------------

def make_fermion_chain(n, seed):
    from dmrgpy import fermionchain
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    mh = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        mh[i, i + 1] = 1.0
        mh[i + 1, i] = 1.0
    for i in range(n):
        mh[i, i] = 1j * (-1) ** i
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + mh[i, j] * fc.Cdag[i] * fc.C[j]
    for i in range(n - 1):
        h = h + (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    fc.set_hamiltonian(h)
    return fc, h


def make_spin_chain(n, seed):
    from dmrgpy import spinchain
    rng = np.random.RandomState(seed)
    spins = ["S=1/2" for i in range(n)]
    sc = spinchain.Spin_Chain(spins)
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * sc.Sx[i] * sc.Sx[i + 1]
        h = h + rng.random() * sc.Sy[i] * sc.Sy[i + 1]
        h = h + rng.random() * sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc, h


def run_case_arnolditk(build, n, mode, seed=1, tol=1e-6, **kwargs):
    chain, h = build(n, seed)
    chain.mode = None
    if mode == "DMRG": chain.setup_python()  # no compiled extension needed
    arnolditk.arnoldimode = mode
    reset_counts()
    t0 = time.time()
    try:
        es, wfs = arnolditk.mpsarnoldi(chain, h, mode="GS", delta=tol, **kwargs)
        ok, err = True, None
    except Exception as e:
        es, wfs = None, None
        ok, err = False, repr(e)
    dt = time.time() - t0
    return dict(ok=ok, err=err,
                energies=(np.atleast_1d(es).tolist() if ok else None),
                time=dt, nops=get_count(mode))


def run_case_iram(build, n, mode, seed=1, tol=1e-6, **kwargs):
    chain, h = build(n, seed)
    chain.mode = None
    if mode == "DMRG": chain.setup_python()  # no compiled extension needed
    arnolditk.arnoldimode = mode  # arpacktk reuses arnolditk's random_state,
                                   # which reads this same module flag
    reset_counts()
    t0 = time.time()
    try:
        es, wfs = arpacktk.mpsiram(chain, h, which="SR", tol=tol, **kwargs)
        ok, err = True, None
    except Exception as e:
        es, wfs = None, None
        ok, err = False, repr(e)
    dt = time.time() - t0
    return dict(ok=ok, err=err,
                energies=(np.atleast_1d(es).tolist() if ok else None),
                time=dt, nops=get_count(mode))


def ground_truth(build, n, seed, nwf=1, n_gt=None):
    chain, h = build(n, seed)
    if n_gt is None:
        n_gt = max(nwf, 6)
    es = chain.get_excited(mode="ED", n=n_gt)
    es = np.array(es)
    order = np.argsort(es.real)
    return es[order].tolist()


CASES = []

def add_case(name, build, n, nwf=1, seed=1, **kwargs):
    CASES.append(dict(name=name, build=build, n=n, nwf=nwf, seed=seed, kwargs=kwargs))

add_case("nonhermitian_fermion_n4", make_fermion_chain, 4, nwf=1, seed=3)
add_case("nonhermitian_fermion_n6", make_fermion_chain, 6, nwf=1, seed=4)
add_case("nonhermitian_fermion_n8", make_fermion_chain, 8, nwf=1, seed=7)
add_case("hermitian_spin_n6", make_spin_chain, 6, nwf=1, seed=1)
add_case("hermitian_spin_n10", make_spin_chain, 10, nwf=1, seed=2)

# Smaller variants for DMRG mode (itensor_version="python"/pyitensor):
# each Op(x) there is a real MPO*MPS application, much slower than an ED
# sparse matvec, so keep chain lengths modest to keep this runnable
# without a compiled C++ extension.
SMALL_CASES = []
def add_small_case(name, build, n, nwf=1, seed=1, **kwargs):
    SMALL_CASES.append(dict(name=name, build=build, n=n, nwf=nwf, seed=seed, kwargs=kwargs))

add_small_case("nonhermitian_fermion_n4", make_fermion_chain, 4, nwf=1, seed=3)
add_small_case("hermitian_spin_n4", make_spin_chain, 4, nwf=1, seed=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--modes", default="ED,DMRG")
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    modes = args.modes.split(",")
    results = []
    all_cases = {"ED": CASES, "DMRG": SMALL_CASES}
    for mode in modes:
        for case in all_cases.get(mode, CASES):
            gt = ground_truth(case["build"], case["n"], case["seed"], nwf=case["nwf"])
            print(f"=== {case['name']} (n={case['n']}, mode={mode}) === ground truth: {np.round(gt[0],6)}")
            for algo, runner in [("arnolditk", run_case_arnolditk),
                                  ("iram", run_case_iram)]:
                res = runner(case["build"], case["n"], mode, tol=args.tol,
                              seed=case["seed"], **case["kwargs"])
                res.update(name=case["name"], mode=mode, algo=algo, n=case["n"])
                if res["ok"]:
                    es = np.array(res["energies"], dtype=complex)
                    gts = np.array(gt, dtype=complex)
                    err = float(np.min(np.abs(es[0] - gts)))
                    res["energies"] = [[e.real, e.imag] for e in es]
                else:
                    err = None
                res["max_energy_error"] = err
                results.append(res)
                print(f"  algo={algo:9s} mode={mode:5s} ok={res['ok']!s:5s} "
                      f"nops={res['nops']:4d} time={res['time']:7.3f}s "
                      f"err={err} "
                      f"energy={res['energies'][0] if res['ok'] else res['err']}")

    print("\n=== Summary: Op(x) call count, arnolditk vs iram ===")
    for mode in modes:
        for case in all_cases.get(mode, CASES):
            a = [r for r in results if r["name"] == case["name"]
                 and r["mode"] == mode and r["algo"] == "arnolditk"][0]
            b = [r for r in results if r["name"] == case["name"]
                 and r["mode"] == mode and r["algo"] == "iram"][0]
            if a["ok"] and b["ok"] and a["nops"] > 0:
                ratio = b["nops"] / a["nops"]
                print(f"  {mode:5s} {case['name']:28s} "
                      f"arnolditk={a['nops']:4d}  iram={b['nops']:4d}  "
                      f"(iram/arnolditk={ratio:.2f}x)")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.json}")

if __name__ == "__main__":
    main()
