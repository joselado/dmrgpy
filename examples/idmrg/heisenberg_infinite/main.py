# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt  # plotting
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#################################################
### Infinite, translationally-invariant chain ###
#################################################
# Infinite_Spin_Chain is defined by a single n_uc-site unit cell (here
# n_uc=1, a uniform chain) instead of a fixed total length -- couplings
# are written with the SxC/SxR (central-cell/next-cell) operator
# convention, see infinitechain.py's own docstring.
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg"

h = (ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0])
ic.set_hamiltonian(h)

ic.maxm = 30      # wavefunction bond dimension cap
ic.maxiter = 40   # iDMRG macro-iterations (growth steps)
ic.etol = 1e-9    # energy-density convergence tolerance
ic.verbose = True  # print the energy density after every macro-iteration

density = ic.gs_energy()

exact = 0.25 - np.log(2)  # Bethe-ansatz ground-state energy density
print()
print("iDMRG converged:", ic.converged)
print("iDMRG energy density:      ", density)
print("exact (Bethe ansatz) value:", exact)
print("absolute error:            ", abs(density - exact))

# Static correlators of the converged infinite chain (transfer-matrix
# formalism, see idmrg.py's onsite_expectation/two_point_correlator).
print()
print("<Sz> (expect ~0 by symmetry):", ic.vev("Sz", 0))
rs = range(1, 8)
correlators = []
for r in rs:
    c = ic.correlator("Sz", 0, "Sz", r)
    correlators.append(c)
    print("<Sz(0)Sz({})> =".format(r), c)

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(list(rs), correlators, "o-", color="tab:blue")
ax.axhline(0, color="gray", lw=0.8, ls=":")
ax.set_xlabel("distance $r$")
ax.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax.set_title("Uniform Heisenberg chain (iDMRG): static Sz-Sz correlator")
fig.tight_layout()
fig.savefig("heisenberg_infinite_correlator.png", dpi=150)
print("\nsaved plot to heisenberg_infinite_correlator.png")
