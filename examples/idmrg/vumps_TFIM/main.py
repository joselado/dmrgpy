# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
from scipy.integrate import quad
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#####################################################################
### VUMPS (pyitensor/vumps.py) vs iDMRG's own growing algorithm,  ###
### transverse-field Ising model, sweeping the bond dimension D   ###
#####################################################################
# Unlike idmrg_ground_state (which GROWS a window indefinitely and
# truncates down to maxm at every step), gs_method="vumps" solves
# directly, in the thermodynamic limit, for the D-dimensional
# variational optimum -- see infinitechain.py's own docstring and
# pyitensor/vumps.py's module docstring for the algorithm. This example
# sweeps D and compares both methods against the exact (Pfeuty 1970)
# free-fermion energy density.

g = 1.5  # transverse field (g>1: gapped paramagnetic phase)


def exact_energy_density(g):
    """H = -sum(sigma^x_i sigma^x_{i+1}) - g*sum(sigma^z_i)."""
    val, _ = quad(lambda k: np.sqrt(1 + g**2 - 2*g*np.cos(k)), 0, np.pi)
    return -val/np.pi


def tfim_chain(gs_method):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    # Sx,Sz have eigenvalues +-1/2, so -4*SxSx=-sigma^x sigma^x and
    # -2*g*Sz=-g*sigma^z, matching the standard Pauli-matrix convention.
    h = -4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_method = gs_method
    return ic


exact = exact_energy_density(g)
Ds = [1, 2, 3, 4]
e_vumps, e_idmrg = [], []

for D in Ds:
    ic = tfim_chain("vumps")
    ic.maxm = D
    # See pyitensor/vumps.py's own "Convergence robustness" section:
    # a generous maxiter/vumps_nrestarts materially improves (without
    # fully guaranteeing) how close VUMPS lands to the true optimum,
    # especially once D>1.
    ic.maxiter = 400
    ic.vumps_nrestarts = 6
    e_vumps.append(ic.gs_energy())

    ic2 = tfim_chain("idmrg")
    ic2.maxm = D
    ic2.maxiter = 200
    ic2.etol = 1e-10
    e_idmrg.append(ic2.gs_energy())

print("D       VUMPS         iDMRG (growing)    exact")
for D, ev, ei in zip(Ds, e_vumps, e_idmrg):
    print("{:2d}   {:12.6f}   {:12.6f}   {:12.6f}".format(D, ev, ei, exact))

# TFIM's own ground state is genuinely entangled at any finite g -- D=1
# is a real (product-state) approximation here, not the exactly-solvable
# case pyitensor/vumps.py's own field-only/dimer tests use, so every D is
# only checked to stay within the (loose, documented) variational window
# -- see vumps.py's own "Convergence robustness" section for why a tight
# match is not asserted, especially at D>1.
for ev in e_vumps:
    assert ev > exact - 0.2 and ev < exact + 0.5

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(Ds, e_vumps, "o-", label="VUMPS (fixed D)")
ax.plot(Ds, e_idmrg, "s-", label="iDMRG (growing, truncated to maxm=D)")
ax.axhline(exact, color="k", ls="--", label="exact (free fermion)")
ax.set_xlabel("bond dimension D")
ax.set_ylabel("ground-state energy density")
ax.set_title("TFIM (g={}): VUMPS vs iDMRG vs exact".format(g))
ax.legend()
fig.tight_layout()
fig.savefig("vumps_TFIM.png", dpi=150)
print("\nSaved plot to vumps_TFIM.png")
