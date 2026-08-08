# Infinite DMRG (iDMRG): pure-Python (pyitensor) backend vs ITensor v3
# (mpscpp3), on the uniform S=1/2 Heisenberg chain, for both a 1-site and
# a 2-site unit cell.
#
# mpscpp3/chain_session.h's Chain::idmrg_ground_state is a line-by-line
# C++ port of pyitensor/idmrg.py's own growing algorithm -- this checks
# both backends land on the same energy density at matched (maxm,
# maxiter, etol), and that a 2-site unit cell reproduces the same
# physical uniform chain as a 1-site unit cell (translational invariance:
# the model doesn't actually care how the unit cell is chosen).
#
# infinitechain.py's public API dispatches itensor_version=3 straight to
# Chain::idmrg_ground_state (energy density only -- no static-correlator
# support on that backend yet, see gs_energy's own docstring); vev()/
# correlator() still require itensor_version="python" regardless of
# which backend gs_energy() itself used.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import infinitechain

MAXM, CUTOFF, MAXITER, ETOL, NITER = 30, 1e-12, 60, 1e-9, 30
EXACT = 0.25 - np.log(2)  # Bethe-ansatz Heisenberg S=1/2 energy density


def build_chain(n_uc, itensor_version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc, itensor_version=itensor_version)
    ic.gs_method = "idmrg"
    if n_uc == 1:
        h = ic.SxC[0] * ic.SxR[0] + ic.SyC[0] * ic.SyR[0] + ic.SzC[0] * ic.SzR[0]
    else:  # n_uc == 2
        h = (ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1] + ic.SzC[0] * ic.SzC[1]
             + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0] + ic.SzC[1] * ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.maxm, ic.cutoff, ic.maxiter, ic.etol, ic.niter = MAXM, CUTOFF, MAXITER, ETOL, NITER
    return ic


results = {}
for n_uc in (1, 2):
    for version in ("python", 3):
        ic = build_chain(n_uc, version)
        t0 = time.time()
        density = ic.gs_energy()
        dt = time.time() - t0
        results[(n_uc, version)] = (density, dt)
        print("n_uc={} itensor_version={:6} density={: .10f}  |err|={:.2e}  time={:6.2f}s".format(
            n_uc, str(version), density, abs(density - EXACT), dt))

print()
for n_uc in (1, 2):
    d_py, t_py = results[(n_uc, "python")]
    d_v3, t_v3 = results[(n_uc, 3)]
    print("n_uc={}: v3 is {:.2f}x python's time; |density diff (python vs v3)|={:.2e}".format(
        n_uc, t_v3 / t_py, abs(d_py - d_v3)))
    assert abs(d_py - d_v3) < 1e-6, "python and v3 iDMRG energy densities disagree"

d_n1 = results[(1, "python")][0]
d_n2 = results[(2, "python")][0]
print("|density diff (n_uc=1 vs n_uc=2, python)|={:.2e}".format(abs(d_n1 - d_n2)))
assert abs(d_n1 - d_n2) < 1e-4, "n_uc=1 and n_uc=2 energy densities disagree"

print()
print("iDMRG python-VS-v3 backend regression test PASSED")

n_ucs = (1, 2)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

for version, marker, label in (("python", "o-", "python"), (3, "s--", "v3")):
    dens = [results[(n_uc, version)][0] for n_uc in n_ucs]
    ax1.plot(n_ucs, dens, marker, label=label)
ax1.axhline(EXACT, color="k", ls=":", label="exact (Bethe ansatz)")
ax1.set_xlabel("unit-cell size $n_{uc}$")
ax1.set_ylabel("ground-state energy density")
ax1.set_xticks(n_ucs)
ax1.set_title("energy density")
ax1.legend()

for version, marker, label in (("python", "o-", "python"), (3, "s--", "v3")):
    times = [results[(n_uc, version)][1] for n_uc in n_ucs]
    ax2.plot(n_ucs, times, marker, label=label)
ax2.set_xlabel("unit-cell size $n_{uc}$")
ax2.set_ylabel("wall time (s)")
ax2.set_xticks(n_ucs)
ax2.set_title("wall time")
ax2.legend()

fig.suptitle("iDMRG: python (pyitensor) vs ITensor v3 (mpscpp3)")
fig.tight_layout()
fig.savefig("heisenberg_infinite_python_VS_v3.png", dpi=150)
print("saved plot to heisenberg_infinite_python_VS_v3.png")
