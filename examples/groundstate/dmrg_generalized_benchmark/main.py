"""Head-to-head benchmark: the self-consistent-Lagrange-multiplier
generalized-eigenvalue DMRG (Many_Body_Chain.gs_energy_generalized())
across its three implementations -- the pure-Python pyitensor engine
(pyitensor/dmrg.py's dmrg_generalized()), its ITensor v3 C++ port
(Chain::gs_energy_generalized, mpscpp3/chain_session.h) and the live
Julia/ITensors.jl one (mpsjulialive/generalized.jl's
get_gs_generalized) -- against the
pre-existing ARPACK-mode-2 route (dmrgpy.algebra.arpacktk.
mpsiram_generalized, OP=inv(M)*A via the approximate correction-vector
self.applyinverse), all solving the same generalized eigenproblem
H|psi>=lambda*A|psi> for a Hermitian positive-definite metric A.

dmrg_generalized() builds an actual materialized MPO H-lambda*A each
outer sweep and runs ordinary two-site DMRG against it, needing no
approximate operator inverse at all; mpsiram_generalized() instead needs
one iterative correction-vector solve (self.applyinverse) per Krylov
step, which is only ever approximate in DMRG mode (there is no exact MPO
inverse -- see arpacktk.py's own module docstring) -- so besides wall
time, this benchmark is also a direct check of how much accuracy that
approximation costs relative to a method with no such step. The three
dmrg_generalized() implementations should agree with each other (and
with ED) to near machine precision -- they run the identical algorithm,
line-for-line, just against a hand-rolled two-site sweep in pure Python
vs the real, compiled ITensor v3 library vs real ITensors.jl -- so the
interesting comparison between them is wall time, not accuracy.

Ground truth is scipy.linalg.eigh's exact generalized Hermitian-definite
eigensolver, built from the ED sparse matrices of H and A (same
convention as tests/test_dmrg_generalized.py and
tests/test_arpacktk_iram.py).

The v3 rows are skipped automatically if the compiled extension isn't
available (see cppext.available(3)) -- run `python install.py
--itensor-version=3` first to include them; likewise the julia_live rows
need a Julia toolchain (`python install_julia.py`).

Ends with a two-panel plot (wall time and accuracy vs chain length) over
every available solver.

Reading the julia_live timings: Julia JIT-compiles on first call, so the
*first* case's wall time is dominated by compilation (~46s vs ~0.1-0.2s
for every later case in a measured run) and says nothing about the
solver. Compare the second and third cases instead -- warm, this backend
sits between the other two at n=8 (v3 0.16s, julia 0.20s, pyitensor
0.50s), with all three at machine-precision accuracy.

Usage: python main.py [--json out.json]
"""
import sys, os, time, json, argparse

sys.path.append(os.getcwd() + '/../../../src')

import numpy as np
import scipy.linalg as sla
import matplotlib.pyplot as plt

from dmrgpy import cppext, fermionchain
from dmrgpy.algebra import arnolditk, arpacktk


# ---------------------------------------------------------------------
# Test problem: an interacting fermion chain Hamiltonian H (random
# hopping + nearest-neighbor interaction) and a diagonal, strictly
# positive metric A = 1 + 0.5*sum(N_i) -- identical construction to
# test_arpacktk_iram.py's _generalized_fermion_problem and
# tests/test_dmrg_generalized.py's own copy, so results here are
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

def run_dmrg_generalized(n, seed, maxm, nsweeps, itensor_version, cutoff=1e-12):
    fc, h, m = generalized_fermion_problem(n, seed)
    if itensor_version == "python":
        fc.setup_python()
    elif itensor_version == "julia_live":
        fc.setup_julia()
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

def julia_available():
    """Whether the live Julia backend can actually be used here (the
    import is what boots the Julia session), mirroring
    cppext.available(3) for the compiled C++ one."""
    try:
        from dmrgpy.mpsjulialive import juliasession  # noqa: F401
        return True
    except Exception:
        return False


V3_AVAILABLE = cppext.available(3)
JULIA_AVAILABLE = julia_available()
ALGOS = ["dmrg_generalized_python"] \
        + (["dmrg_generalized_v3"] if V3_AVAILABLE else []) \
        + (["dmrg_generalized_julia"] if JULIA_AVAILABLE else []) \
        + ["arpack_mode2"]


def run_algo(algo, case, tol):
    if algo == "dmrg_generalized_python":
        return run_dmrg_generalized(case["n"], case["seed"], case["maxm"],
                case["nsweeps"], "python")
    if algo == "dmrg_generalized_v3":
        return run_dmrg_generalized(case["n"], case["seed"], case["maxm"],
                case["nsweeps"], 3)
    if algo == "dmrg_generalized_julia":
        return run_dmrg_generalized(case["n"], case["seed"], case["maxm"],
                case["nsweeps"], "julia_live")
    if algo == "arpack_mode2":
        return run_arpack_mode2(case["n"], case["seed"], tol=tol)
    raise ValueError(algo)


def plot(results):
    """Wall time and accuracy vs chain length, one line per solver."""
    fig, (ax_t, ax_e) = plt.subplots(1, 2, figsize=(11, 4.2))
    for algo in ALGOS:
        rows = [r for r in results if r["algo"] == algo and r["ok"]]
        if not rows:
            continue
        ns = [r["n"] for r in rows]
        ax_t.plot(ns, [r["time"] for r in rows], "o-", label=algo)
        # errors can hit exactly 0 at machine precision; floor them so the
        # log axis still shows the point
        ax_e.plot(ns, [max(r["error"], 1e-16) for r in rows], "o-", label=algo)
    ax_t.set_xlabel("chain length n")
    ax_t.set_ylabel("wall time [s]")
    ax_t.set_yscale("log")
    ax_t.set_title("Generalized eigenproblem: wall time")
    ax_e.set_xlabel("chain length n")
    ax_e.set_ylabel(r"$|\lambda-\lambda_\mathrm{exact}|$")
    ax_e.set_yscale("log")
    ax_e.set_title("Accuracy vs scipy.linalg.eigh")
    for ax in (ax_t, ax_e):
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--tol", type=float, default=1e-8,
            help="ARPACK mode-2 convergence tolerance")
    args = ap.parse_args()

    if not V3_AVAILABLE:
        print("(ITensor v3 extension not compiled -- skipping "
              "dmrg_generalized_v3 rows; run `python install.py "
              "--itensor-version=3` to include them)\n")
    if not JULIA_AVAILABLE:
        print("(no Julia toolchain -- skipping dmrg_generalized_julia "
              "rows; run `python install_julia.py` to include them)\n")

    results = []
    for case in CASES:
        fc, h, m = generalized_fermion_problem(case["n"], case["seed"])
        gt = ground_truth(fc, h, m)
        print(f"=== {case['name']} (n={case['n']}) === ground truth lambda0 = {gt:.8f}")

        for algo in ALGOS:
            res = run_algo(algo, case, args.tol)
            err = abs(res["energy"] - gt) if res["ok"] else None
            res.update(name=case["name"], n=case["n"], algo=algo, ground_truth=gt,
                    error=err)
            results.append(res)
            e_str = f"{res['energy']:.8f}" if res["ok"] else res["err"]
            err_str = f"{err:.2e}" if err is not None else "n/a"
            print(f"  algo={algo:24s} ok={res['ok']!s:5s} "
                  f"time={res['time']:7.3f}s  error={err_str:>10s}  energy={e_str}")

    print("\n=== Summary: wall time relative to dmrg_generalized_python ===")
    for case in CASES:
        rows = {r["algo"]: r for r in results if r["name"] == case["name"]}
        base = rows["dmrg_generalized_python"]
        print(f"  {case['name']}:")
        for algo in ALGOS:
            r = rows[algo]
            if not (r["ok"] and base["ok"]):
                continue
            ratio = r["time"] / base["time"] if base["time"] > 0 else float("inf")
            print(f"    {algo:24s} {r['time']:7.3f}s (err {r['error']:.2e})"
                  f"   -> {ratio:.2f}x the pure-Python dmrg_generalized time")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.json}")

    plot(results)


if __name__ == "__main__":
    main()
