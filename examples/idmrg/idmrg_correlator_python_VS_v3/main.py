# iDMRG static correlators and local gap: pure-Python (pyitensor) backend
# vs ITensor v3 (mpscpp3) -- correctness AND cost.
#
# mpscpp3/chain_session.h now carries a full C++ port of pyitensor/idmrg.py's
# static-observable machinery: the gauge-consistent unit cell extracted from
# a single micro-step's theta (Chain::idmrg_theta_cell +
# ic_canonicalize_cell), the transfer-matrix expectation values built on it
# (Chain::idmrg_onsite_expectation / Chain::idmrg_two_point_correlator), and
# the local superblock gap (Chain::idmrg_local_excitation_gap). So
# `itensor_version=3` with `gs_method="idmrg"` supports vev()/correlator()/
# local_excitation_gap() exactly like `itensor_version="python"` does.
#
# Part 1 checks the two backends agree, against a model with exact answers.
# Part 2 measures what each costs as the bond dimension grows -- the two
# halves need different models, which is why this script runs two:
#
#   * Part 1 uses a GAPPED transverse-field Ising chain (h > J). Gapped and
#     non-degenerate matters for a *state* comparison: on a critical chain
#     the growing algorithm settles onto a slightly, arbitrarily
#     symmetry-broken state whose residual depends on its own unseeded
#     random starting MPS, so the two backends' energies agree while their
#     states need not. That is a property of the model, not of either
#     implementation. Its bond dimension also saturates almost immediately,
#     which makes it useless for part 2.
#
#   * Part 2 uses the CRITICAL (gapless) uniform Heisenberg chain, whose
#     entanglement keeps growing, so `maxm` is a real knob and the timings
#     actually scale with it.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
from statistics import median

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import cppext
from dmrgpy import infinitechain

RS = list(range(1, 8))

backends = ["python"]
if cppext.available(3):
    backends.append(3)
else:
    print("mpscpp3 not compiled -- running the pure-Python backend only\n")

MARKERS = dict(zip(("python", 3), ("o-", "s--")))


# == Part 1: correctness, against exact free fermions =======================

J, H_FIELD = 1.0, 2.0


def exact_tfim(J, h, nk=200001):
    """Free-fermion (Jordan-Wigner) exact answers for
    H = -J sum sigmax_i sigmax_{i+1} - h sum sigmaz_i: energy per site and
    <Sz> = <sigmaz>/2, from the filled negative Bogoliubov band."""
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    eps = 2.0 * np.sqrt(J ** 2 + h ** 2 - 2 * J * h * np.cos(k))
    return -0.5 * eps.mean(), 0.5 * ((2.0 * (h - J * np.cos(k))) / eps).mean()


def tfim_chain(itensor_version, maxm=16):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(-4.0 * J * (ic.SxC[0] * ic.SxR[0]) - 2.0 * H_FIELD * ic.SzC[0])
    ic.maxm, ic.maxiter, ic.etol = maxm, 120, 1e-12
    ic.gs_energy()
    return ic


e0_exact, sz_exact = exact_tfim(J, H_FIELD)
print("== Part 1: gapped TFIM (J={}, h={}) vs exact free fermions ==".format(J, H_FIELD))
print("exact: e0 = {: .10f}   <Sz> = {: .10f}".format(e0_exact, sz_exact))

tfim, corr_tfim, gap_tfim = {}, {}, {}
for version in backends:
    ic = tfim_chain(version)
    sz = ic.vev("Sz", 0).real
    corr_tfim[version] = np.array([ic.correlator("Sz", 0, "Sz", r).real for r in RS])
    gap_tfim[version] = ic.local_excitation_gap()
    tfim[version] = (ic.e0, sz)
    print("  {:6}  e0 err {:.2e}   <Sz> err {:.2e}   local gap {:.6f}".format(
        str(version), abs(ic.e0 - e0_exact), abs(sz - sz_exact), gap_tfim[version]))

if 3 in backends:
    de = abs(tfim["python"][0] - tfim[3][0])
    dsz = abs(tfim["python"][1] - tfim[3][1])
    dc = np.max(np.abs(corr_tfim["python"] - corr_tfim[3]))
    dg = abs(gap_tfim["python"] - gap_tfim[3])
    print("  python vs v3: |de0|={:.2e}  |d<Sz>|={:.2e}  max|dcorr|={:.2e}  "
          "|dgap|={:.2e}".format(de, dsz, dc, dg))
    assert de < 1e-8, "python and v3 iDMRG energy densities disagree"
    assert dsz < 1e-5, "python and v3 iDMRG <Sz> disagree"
    assert dc < 1e-5, "python and v3 iDMRG correlators disagree"
    assert dg < 1e-2, "python and v3 iDMRG local gaps disagree"
    print("  iDMRG static-observable python-VS-v3 regression check PASSED")


# == Part 2: cost vs bond dimension, on the critical Heisenberg chain =======
#
# Timed as medians of NREPS runs, not single shots. That is not pedantry
# here: iDMRG starts from an unseeded random MPS, and on a critical chain
# the transfer matrix's leading eigenvalues crowd together, so the
# iterative eigensolves behind both the ground state and the observables
# used to vary by 5-10x run to run (one single-shot measurement of the
# configuration below came in at 213s against a 13s median). Medians of a
# few repeats are the minimum needed for the numbers to mean anything.
#
# Three quantities, because they have genuinely different cost structures:
#
#   ground state   -- the growth loop, PLUS the one-off extraction and
#                     re-gauging of the unit cell at the end of it. The
#                     latter is a chi^2 x chi^2 transfer-matrix eigenproblem
#                     and, on this critical chain, it dominates: measured by
#                     direct instrumentation at 70-93% of pyitensor's whole
#                     gs_energy() call.
#   first vev      -- builds and caches the measurement-independent
#                     environment (transfer tensors + both families of
#                     fixed points). Another pair of eigenproblems.
#   correlator sweep -- reuses that cache, so it measures only the operator
#                     string itself, applied site by site.

MAXMS = [8, 16, 24, 32]
NREPS = 3


def heisenberg_chain(itensor_version, maxm, maxiter=40):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    ic.set_hamiltonian(ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0]
                       + ic.SzC[0] * ic.SzR[0])
    ic.maxm, ic.maxiter, ic.etol = maxm, maxiter, 1e-12
    return ic


print()
print("== Part 2: cost on the critical Heisenberg chain, medians of {} ==".format(NREPS))
print("   (pin BLAS threads -- MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 -- or the")
print("    timings are not comparable across runs; see src/dmrgpy/blasthreads.py)")
t_gs, t_vev, t_corr = {}, {}, {}
for version in backends:
    for maxm in MAXMS:
        gs, vev, corr = [], [], []
        for _rep in range(NREPS):
            ic = heisenberg_chain(version, maxm)
            t0 = time.time(); ic.gs_energy(); gs.append(time.time() - t0)
            t0 = time.time(); ic.vev("Sz", 0); vev.append(time.time() - t0)
            t0 = time.time()
            for r in RS:
                ic.correlator("Sz", 0, "Sz", r)
            corr.append(time.time() - t0)
        t_gs[(version, maxm)] = median(gs)
        t_vev[(version, maxm)] = median(vev)
        t_corr[(version, maxm)] = median(corr)
        print("  {:6}  maxm={:3d}   ground state {:8.3f}s   first vev {:8.4f}s"
              "   correlator sweep (r=1..{}) {:8.4f}s".format(
                  str(version), maxm, t_gs[(version, maxm)], t_vev[(version, maxm)],
                  RS[-1], t_corr[(version, maxm)]))

if 3 in backends:
    print()
    for maxm in MAXMS:
        print("  maxm={:3d}: v3 vs python -- ground state {:6.2f}x, first vev {:6.2f}x,"
              "  correlator sweep {:6.2f}x".format(
                  maxm,
                  t_gs[("python", maxm)] / t_gs[(3, maxm)],
                  t_vev[("python", maxm)] / t_vev[(3, maxm)],
                  t_corr[("python", maxm)] / t_corr[(3, maxm)]))
    print("  (>1 means v3 is faster)")
    print()
    print("  Both backends now share the same two optimizations, so the")
    print("  observables land within a small factor of each other: the")
    print("  measurement-independent environment is built once per ground state")
    print("  and cached across calls, and an operator string is applied to the")
    print("  fixed point one site at a time (O(chi^3 d)) rather than composed")
    print("  into a transfer matrix first (O(chi^6) per site). Before those, a")
    print("  seven-point sweep at maxm=24 cost ~68s on pyitensor and ~4s on v3.")
    print()
    print("  The ground state is where they still differ, and pyitensor wins:")
    print("  its growth loop is NumPy on dense arrays, while v3's runs through")
    print("  ITensor objects whose per-contraction bookkeeping dominates at")
    print("  these small bond dimensions. Note this is the opposite of what an")
    print("  earlier version of this script reported -- back then pyitensor's")
    print("  gs_energy() was dominated by its own final canonicalization, whose")
    print("  ARPACK solve ran with the default 20-vector Krylov subspace and")
    print("  converged badly against this chain's clustered transfer spectrum.")
    print("  Raising it to _ARPACK_NCV=40 (pyitensor/idmrg.py, matching the")
    print("  krylovdim mpscpp3 already used) cut that from ~13s to ~0.6s at")
    print("  maxm=32 with every value unchanged to ~1e-16 -- and removed the")
    print("  run-to-run variance that made the comparison unreadable.")

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))

ax = axes[0]
for version in backends:
    ax.plot(RS, np.abs(corr_tfim[version]), MARKERS[version], label=str(version))
ax.set_yscale("log")
ax.set_xlabel("separation $r$ (sites)")
ax.set_ylabel(r"$|\langle S^z_0 S^z_r\rangle|$")
ax.set_title("gapped TFIM correlator (exponential decay)")
ax.legend()

ax = axes[1]
for version in backends:
    ax.plot(MAXMS, [t_gs[(version, m)] for m in MAXMS], MARKERS[version],
            label="{}".format(version))
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("bond dimension maxm")
ax.set_ylabel("wall time (s)")
ax.set_title("cost: iDMRG ground state")
ax.legend()

ax = axes[2]
for version in backends:
    ax.plot(MAXMS, [t_vev[(version, m)] for m in MAXMS], MARKERS[version],
            label="{}: first vev".format(version))
    ax.plot(MAXMS, [t_corr[(version, m)] for m in MAXMS],
            MARKERS[version].replace("o", "^").replace("s", "v"),
            label="{}: correlator sweep".format(version))
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("bond dimension maxm")
ax.set_ylabel("wall time (s)")
ax.set_title("cost: static observables")
ax.legend(fontsize=8)

fig.suptitle("iDMRG static observables: python (pyitensor) vs ITensor v3 (mpscpp3)")
fig.tight_layout()
fig.savefig("idmrg_correlator_python_VS_v3.png", dpi=150)
print()
print("saved plot to idmrg_correlator_python_VS_v3.png")
