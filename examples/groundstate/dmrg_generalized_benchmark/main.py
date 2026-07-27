"""Head-to-head benchmark: the new self-consistent-Lagrange-multiplier
generalized-eigenvalue DMRG (Many_Body_Chain.gs_energy_generalized(),
pyitensor/dmrg.py's dmrg_generalized() -- see its own docstring for the
algorithm) vs the pre-existing ARPACK-mode-2 route
(dmrgpy.algebra.arpacktk.mpsiram_generalized, OP=inv(M)*A via the
approximate correction-vector self.applyinverse), both solving the same
generalized eigenproblem H|psi>=lambda*A|psi> for a Hermitian
positive-definite metric A, in DMRG mode on the pure-Python (pyitensor)
backend -- the only backend the new solver is implemented on so far (see
CLAUDE.md's "Implement in pyitensor first" scoping).

dmrg_generalized() builds an actual materialized MPO H-lambda*A each
outer sweep and runs ordinary two-site DMRG against it, needing no
approximate operator inverse at all; mpsiram_generalized() instead needs
one iterative correction-vector solve (self.applyinverse) per Krylov
step, which is only ever approximate in DMRG mode (there is no exact MPO
inverse -- see arpacktk.py's own module docstring) -- so besides wall
time, this benchmark is also a direct check of how much accuracy that
approximation costs relative to a method with no such step.

Ground truth is scipy.linalg.eigh's exact generalized Hermitian-definite
eigensolver, built from the ED sparse matrices of H and A (same
convention as tests/test_arpacktk_iram.py and
tests/test_pyitensor_dmrg_generalized.py).

Usage: python main.py [--json out.json]
"""
import sys, os, time, json, argparse

sys.path.append(os.getcwd() + '/../../../src')

import numpy as np
import scipy.linalg as sla

from dmrgpy import fermionchain
from dmrgpy.algebra import arnolditk, arpacktk


# ---------------------------------------------------------------------
# Test problem: an interacting fermion chain Hamiltonian H (random
# hopping + nearest-neighbor interaction) and a diagonal, strictly
# positive metric A = 1 + 0.5*sum(N_i) -- identical construction to
# test_arpacktk_iram.py's _generalized_fermion_problem and
# test_pyitensor_dmrg_generalized.py's own copy, so results here are
# directly comparable to those tests' tolerances.
# ---------------------------------------------------------------------

def generalized_fermion_problem(n, seed):
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    for i in range(n - 1):
        h = h + rng.random() * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    m = 1
    for i in range(n):
        m = m + 0.5 * fc.N[i]
    fc.set_hamiltonian(h)
    return fc, h, m


def ground_truth(fc, h, m):
    Hmat = fc.get_ED_obj().MO2matrix(h)
    Mmat = fc.get_ED_obj().MO2matrix(m)
    Hmat = Hmat.toarray() if hasattr(Hmat, "toarray") else np.array(Hmat)
    Mmat = Mmat.toarray() if hasattr(Mmat, "toarray") else np.array(Mmat)
    return np.sort(sla.eigh(Hmat, Mmat, eigvals_only=True))[0]


# ---------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------

def run_dmrg_generalized(n, seed, maxm, nsweeps, cutoff=1e-12):
    fc, h, m = generalized_fermion_problem(n, seed)
    fc.setup_python()
    fc.maxm, fc.nsweeps, fc.cutoff = maxm, nsweeps, cutoff
    t0 = time.time()
    try:
        lam = fc.gs_energy_generalized(m)
        ok, err = True, None
    except Exception as e:
        lam, ok, err = None, False, repr(e)
    dt = time.time() - t0
    return dict(ok=ok, err=err, energy=lam, time=dt)


def run_arpack_mode2(n, seed, tol=1e-8, **kwargs):
    fc, h, m = generalized_fermion_problem(n, seed)
    fc.setup_python()
    arnolditk.arnoldimode = "DMRG"
    t0 = time.time()
    try:
        es, wfs = arpacktk.mpsiram_generalized(fc, h, m, which="SR", nev=1,
                tol=tol, **kwargs)
        lam, ok, err = es[0].real, True, None
    except Exception as e:
        lam, ok, err = None, False, repr(e)
    dt = time.time() - t0
    return dict(ok=ok, err=err, energy=lam, time=dt)


CASES = [
    dict(name="fermion_n4", n=4, seed=2, maxm=40, nsweeps=12),
    dict(name="fermion_n6", n=6, seed=2, maxm=60, nsweeps=14),
    dict(name="fermion_n8", n=8, seed=3, maxm=80, nsweeps=16),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--tol", type=float, default=1e-8,
            help="ARPACK mode-2 convergence tolerance")
    args = ap.parse_args()

    results = []
    for case in CASES:
        fc, h, m = generalized_fermion_problem(case["n"], case["seed"])
        gt = ground_truth(fc, h, m)
        print(f"=== {case['name']} (n={case['n']}) === ground truth lambda0 = {gt:.8f}")

        r_dmrg = run_dmrg_generalized(case["n"], case["seed"], case["maxm"], case["nsweeps"])
        r_iram = run_arpack_mode2(case["n"], case["seed"], tol=args.tol)

        for algo, res in [("dmrg_generalized", r_dmrg), ("arpack_mode2", r_iram)]:
            err = abs(res["energy"] - gt) if res["ok"] else None
            res.update(name=case["name"], n=case["n"], algo=algo, ground_truth=gt,
                    error=err)
            results.append(res)
            e_str = f"{res['energy']:.8f}" if res["ok"] else res["err"]
            err_str = f"{err:.2e}" if err is not None else "n/a"
            print(f"  algo={algo:17s} ok={res['ok']!s:5s} "
                  f"time={res['time']:7.3f}s  error={err_str:>10s}  energy={e_str}")

    print("\n=== Summary: wall time and accuracy, dmrg_generalized vs arpack_mode2 ===")
    for case in CASES:
        a = [r for r in results if r["name"] == case["name"] and r["algo"] == "dmrg_generalized"][0]
        b = [r for r in results if r["name"] == case["name"] and r["algo"] == "arpack_mode2"][0]
        if a["ok"] and b["ok"]:
            ratio = b["time"] / a["time"] if a["time"] > 0 else float("inf")
            faster = "dmrg_generalized" if ratio >= 1 else "arpack_mode2"
            print(f"  {case['name']:12s} dmrg_generalized: {a['time']:7.3f}s (err {a['error']:.2e})"
                  f"   arpack_mode2: {b['time']:7.3f}s (err {b['error']:.2e})"
                  f"   -> {faster} is {max(ratio, 1/ratio):.2f}x faster")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.json}")


if __name__ == "__main__":
    main()
