"""Head-to-head benchmark: the non-Hermitian generalized-eigenvalue
NH-DMRG solver (Many_Body_Chain.gs_energy_generalized() dispatching to
nhdmrg_generalized() for a non-Hermitian H) across its two
implementations -- pyitensor (pyitensor/nhdmrg.py) and its ITensor v3 C++
port (Chain::nhdmrg_generalized, mpscpp3/chain_session.h) -- against the
pre-existing ARPACK-mode-2 route (dmrgpy.algebra.arpacktk.
mpsiram_generalized, OP=inv(M)*A via the approximate correction-vector
self.applyinverse), all solving H|psi_R>=lambda*A|psi_R> for a possibly
non-Hermitian H and a Hermitian positive-definite metric A.

Unlike the Hermitian generalized-eigenvalue benchmark
(examples/groundstate/dmrg_generalized_benchmark), ARPACK mode 2 here
needs no adaptation at all for a non-Hermitian primary operator: its own
docstring only ever required M (the metric) to be Hermitian
positive-definite for the M-weighted Krylov inner product -- A itself was
always allowed to be general, and hessenberg_ritz's np.linalg.eig (not
eigh) already handles a non-Hermitian/complex spectrum. So this
comparison is a fair like-for-like: both solvers target the eigenvalue
of smallest real part (this codebase's "ground state" convention for
non-Hermitian operators throughout), and nhdmrg_generalized() is a direct
non-Hermitian generalization of the Hermitian dmrg_generalized() using
the exact same self-consistent Lagrange-multiplier trick -- see its own
docstring (complex lambda, biorthogonal Rayleigh quotient
<psi_L|H|psi_R>/<psi_L|A|psi_R> in place of the plain <psi|H|psi>/
<psi|A|psi>).

Ground truth is scipy.linalg.eig's general (non-symmetric) generalized
eigensolver, built from the ED sparse matrices of H and A (same
convention as tests/test_nhdmrg_generalized.py).

The v3 rows are skipped automatically if the compiled extension isn't
available (see cppext.available(3)) -- run `python install.py
--itensor-version=3` first to include them.

Usage: python main.py [--json out.json]
"""
import sys, os, time, json, argparse

sys.path.append(os.getcwd() + '/../../../src')

import numpy as np
import scipy.linalg as sla

from dmrgpy import cppext, fermionchain
from dmrgpy.algebra import arnolditk, arpacktk


# ---------------------------------------------------------------------
# Test problem: a non-Hermitian interacting fermion chain (staggered
# imaginary on-site potential, the model of examples/non_hermitian/
# non_hermitian_chain and tests/test_nh_dmrg.py's own nh_fermion_chain)
# and a diagonal, strictly positive metric A = 1 + 0.5*sum(N_i).
# ---------------------------------------------------------------------

def nh_generalized_fermion_problem(n, seed):
    rng = np.random.RandomState(seed)
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n - 1):
        h = h + rng.random() * (fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i])
    for i in range(n):
        h = h + 1j * (-1) ** i * 0.3 * fc.N[i]
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
    w = sla.eig(Hmat, Mmat, right=False)
    return w[np.argsort(w.real)][0]


# ---------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------

def run_nhdmrg_generalized(n, seed, maxm, nsweeps, itensor_version, cutoff=1e-12):
    fc, h, m = nh_generalized_fermion_problem(n, seed)
    if itensor_version == "python":
        fc.setup_python()
    else:
        fc.setup_cpp(version=itensor_version)
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
    fc, h, m = nh_generalized_fermion_problem(n, seed)
    fc.setup_python()
    arnolditk.arnoldimode = "DMRG"
    t0 = time.time()
    try:
        es, wfs = arpacktk.mpsiram_generalized(fc, h, m, which="SR", nev=1,
                tol=tol, **kwargs)
        lam, ok, err = es[0], True, None
    except Exception as e:
        lam, ok, err = None, False, repr(e)
    dt = time.time() - t0
    return dict(ok=ok, err=err, energy=lam, time=dt)


CASES = [
    dict(name="nh_fermion_n4", n=4, seed=2, maxm=40, nsweeps=14),
    dict(name="nh_fermion_n6", n=6, seed=2, maxm=60, nsweeps=16),
    dict(name="nh_fermion_n8", n=8, seed=3, maxm=80, nsweeps=18),
]

V3_AVAILABLE = cppext.available(3)
ALGOS = ["nhdmrg_generalized_python"] + (["nhdmrg_generalized_v3"] if V3_AVAILABLE else []) \
        + ["arpack_mode2"]


def run_algo(algo, case, tol):
    if algo == "nhdmrg_generalized_python":
        return run_nhdmrg_generalized(case["n"], case["seed"], case["maxm"],
                case["nsweeps"], "python")
    if algo == "nhdmrg_generalized_v3":
        return run_nhdmrg_generalized(case["n"], case["seed"], case["maxm"],
                case["nsweeps"], 3)
    if algo == "arpack_mode2":
        return run_arpack_mode2(case["n"], case["seed"], tol=tol)
    raise ValueError(algo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--tol", type=float, default=1e-8,
            help="ARPACK mode-2 convergence tolerance")
    args = ap.parse_args()

    if not V3_AVAILABLE:
        print("(ITensor v3 extension not compiled -- skipping "
              "nhdmrg_generalized_v3 rows; run `python install.py "
              "--itensor-version=3` to include them)\n")

    results = []
    for case in CASES:
        fc, h, m = nh_generalized_fermion_problem(case["n"], case["seed"])
        gt = ground_truth(fc, h, m)
        print(f"=== {case['name']} (n={case['n']}) === ground truth lambda0 = "
              f"{gt.real:.8f}{gt.imag:+.2e}j")

        for algo in ALGOS:
            res = run_algo(algo, case, args.tol)
            err = abs(res["energy"] - gt) if res["ok"] else None
            res.update(name=case["name"], n=case["n"], algo=algo,
                    ground_truth=[gt.real, gt.imag], error=err)
            results.append(res)
            if res["ok"]:
                e = res["energy"]
                e_str = f"{e.real:.8f}{e.imag:+.2e}j"
            else:
                e_str = res["err"]
            err_str = f"{err:.2e}" if err is not None else "n/a"
            print(f"  algo={algo:26s} ok={res['ok']!s:5s} "
                  f"time={res['time']:7.3f}s  error={err_str:>10s}  energy={e_str}")

    print("\n=== Summary: wall time relative to nhdmrg_generalized_python ===")
    for case in CASES:
        rows = {r["algo"]: r for r in results if r["name"] == case["name"]}
        base = rows["nhdmrg_generalized_python"]
        print(f"  {case['name']}:")
        for algo in ALGOS:
            r = rows[algo]
            if not (r["ok"] and base["ok"]):
                continue
            ratio = r["time"] / base["time"] if base["time"] > 0 else float("inf")
            print(f"    {algo:26s} {r['time']:7.3f}s (err {r['error']:.2e})"
                  f"   -> {ratio:.2f}x the pure-Python nhdmrg_generalized time")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.json}")


if __name__ == "__main__":
    main()
