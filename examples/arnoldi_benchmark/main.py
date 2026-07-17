"""Benchmark harness for the arnolditk (Arnoldi) algorithm.

Instruments the number of expensive Op(x) calls (MPO*MPS application in
DMRG mode / sparse matvec in ED mode) alongside wall time and accuracy
against an ED ground truth, for a handful of representative cases:
Hermitian spin-chain ground state, non-Hermitian fermion-chain ground
state, and multi-excited-state cases (both the simultaneous and the
recursive/deflated variants). Since ED is exact for these small chains,
comparing ED vs DMRG results here is also a correctness cross-check of
the Arnoldi implementation itself, not just a performance benchmark.

Usage: python main.py [--json out.json] [--modes ED,DMRG]
"""
import sys, os, time, json, argparse

sys.path.append(os.getcwd()+'/../../src')

import numpy as np

from dmrgpy.multioperatortk.staticoperator import StaticOperator
from dmrgpy.edtk.edchain import EDOperator, State
from dmrgpy.mps import MPS
from dmrgpy.algebra import arnolditk

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
# Test cases
# ---------------------------------------------------------------------

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


def run_case(name, build, n, mode, itensor_version=2, seed=1, **arnoldi_kwargs):
    """Run one Arnoldi case in the given mode, return a result dict."""
    chain, h = build(n, seed)
    chain.itensor_version = itensor_version
    chain.mode = None  # let arnolditk's own arnoldimode flag decide backend
    arnolditk.arnoldimode = mode

    reset_counts()
    t0 = time.time()
    try:
        es, wfs = arnolditk.mpsarnoldi(chain, h, mode="GS", **arnoldi_kwargs)
        ok = True
        err = None
    except Exception as e:
        es, wfs = None, None
        ok = False
        err = repr(e)
    dt = time.time() - t0
    nops = get_count(mode)

    return dict(name=name, mode=mode, n=n, ok=ok, err=err,
                energies=(np.atleast_1d(es).tolist() if ok else None),
                time=dt, nops=nops)


def ground_truth(build, n, seed, nwf=1, n_gt=None):
    """Return a handful of low-lying ED eigenvalues (more than nwf) so
    that matching against Arnoldi's output can be done by nearest-neighbor
    distance rather than assumed ordering -- needed because mode="GS"
    only targets minimum Re(E), and near-degenerate/conjugate-pair
    eigenvalues can be picked up in either branch depending on the random
    Krylov start."""
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

add_case("hermitian_spin_n6", make_spin_chain, 6, nwf=1, seed=1)
add_case("hermitian_spin_n10", make_spin_chain, 10, nwf=1, seed=2)
add_case("nonhermitian_fermion_n4", make_fermion_chain, 4, nwf=1, seed=3)
add_case("nonhermitian_fermion_n6", make_fermion_chain, 6, nwf=1, seed=4)
add_case("multi_excited_spin_n6", make_spin_chain, 6, nwf=3, seed=5,
          recursive_arnoldi=False)
add_case("multi_excited_spin_n6_recursive", make_spin_chain, 6, nwf=3, seed=5,
          recursive_arnoldi=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--modes", default="ED,DMRG")
    args = ap.parse_args()

    modes = args.modes.split(",")
    results = []
    for case in CASES:
        gt = ground_truth(case["build"], case["n"], case["seed"], nwf=case["nwf"])
        print(f"=== {case['name']} (n={case['n']}, nwf={case['nwf']}) === ground truth: {np.round(gt,6)}")
        for mode in modes:
            kwargs = dict(case["kwargs"])
            if case["nwf"] > 1:
                kwargs["nwf"] = case["nwf"]
            res = run_case(case["name"], case["build"], case["n"], mode,
                            seed=case["seed"], **kwargs)
            res["ground_truth"] = [[complex(x).real, complex(x).imag] for x in gt]
            if res["ok"]:
                es = np.array(res["energies"], dtype=complex)
                gts = np.array(gt, dtype=complex)
                # nearest-neighbor distance, robust to conjugate-pair
                # ambiguity in near-degenerate non-Hermitian spectra
                err = float(np.max([np.min(np.abs(e - gts)) for e in es]))
                res["energies"] = [[e.real, e.imag] for e in es]
                if len(es) > 1:
                    pd = [abs(es[i] - es[j]) for i in range(len(es))
                          for j in range(i + 1, len(es))]
                    res["min_pairwise_distinct"] = float(np.min(pd))
                else:
                    res["min_pairwise_distinct"] = None
            else:
                err = None
                res["min_pairwise_distinct"] = None
            res["max_energy_error"] = err
            results.append(res)
            print(f"  mode={mode:5s} ok={res['ok']!s:5s} "
                  f"nops={res['nops']:4d} time={res['time']:7.3f}s "
                  f"err={err} min_pairwise_distinct={res['min_pairwise_distinct']} "
                  f"energies={np.round(res['energies'] if not res['ok'] else np.array(res['energies']),6) if res['ok'] else res['err']}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.json}")

if __name__ == "__main__":
    main()
