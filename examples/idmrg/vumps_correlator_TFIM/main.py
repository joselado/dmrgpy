# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

#####################################################################
### Static correlators from a VUMPS ground state (gs_method=      ###
### "vumps"), compared against iDMRG's own growing algorithm      ###
### (gs_method="idmrg") -- transverse-field Ising model.          ###
#####################################################################
# pyitensor/vumps.py's onsite_expectation/two_point_correlator compute
# <O> and <O_i O_j> directly from the converged mixed-gauge {AC, AR} of
# a VUMPSResult (Vanderstraeten, Haegeman, Verstraete,
# "Tangent-space methods for uniform matrix product states",
# arXiv:1810.07006, Eq.(34)/(37)-(39)) -- no eigenproblem needed, unlike
# pyitensor/idmrg.py's own dominant-right-fixed-point-based correlators,
# since AL/AR are already exactly canonical by construction of the
# VUMPS fixed point. Both are wired into Infinite_Many_Body_Chain.vev/
# correlator, dispatching on ic.gs_method -- this example runs both at
# the same target bond dimension D and compares.

g = 1.5  # transverse field (g>1: gapped paramagnetic phase)
D = 3


def tfim_chain(gs_method):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    # Sx,Sz have eigenvalues +-1/2, so -4*SxSx=-sigma^x sigma^x and
    # -2*g*Sz=-g*sigma^z, matching the standard Pauli-matrix TFIM
    # convention H=-sigma^x sigma^x - g*sigma^z.
    h = -4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0]
    ic.set_hamiltonian(h)
    ic.gs_method = gs_method
    ic.maxm = D
    return ic


ic_vumps = tfim_chain("vumps")
ic_vumps.maxiter = 800
ic_vumps.vumps_nrestarts = 6
ic_vumps.gs_energy()

ic_idmrg = tfim_chain("idmrg")
ic_idmrg.maxiter = 200
ic_idmrg.etol = 1e-10
ic_idmrg.gs_energy()

print("VUMPS  e0={:.8f} converged={}".format(ic_vumps.e0, ic_vumps.converged))
print("iDMRG  e0={:.8f} converged={}".format(ic_idmrg.e0, ic_idmrg.converged))

print()
print("<Sz>: VUMPS={:.6f}  iDMRG={:.6f}".format(
    ic_vumps.vev("Sz", 0).real, ic_idmrg.vev("Sz", 0).real))

rs = list(range(0, 8))
corr_vumps = [ic_vumps.correlator("Sz", 0, "Sz", r).real for r in rs]
corr_idmrg = [ic_idmrg.correlator("Sz", 0, "Sz", r).real for r in rs]

print()
print(" r    <SzSz> (VUMPS)   <SzSz> (iDMRG)")
for r, cv, ci in zip(rs, corr_vumps, corr_idmrg):
    print("{:2d}   {:14.6f}   {:14.6f}".format(r, cv, ci))

# Both methods target the SAME D-dimensional variational state (up to
# each solver's own numerical search), so their correlators should agree
# -- a generous tolerance since VUMPS's own D>1 convergence robustness is
# a documented, open limitation (pyitensor/vumps.py's own "Convergence
# robustness" docstring section) and the two solvers need not land on
# bit-identical states even at the same D.
for cv, ci in zip(corr_vumps, corr_idmrg):
    assert abs(cv - ci) < 0.1

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(rs, corr_vumps, "o-", label="VUMPS (gs_method=\"vumps\")")
ax.plot(rs, corr_idmrg, "s--", label="iDMRG (gs_method=\"idmrg\")")
ax.axhline(0, color="gray", lw=0.8, ls=":")
ax.set_xlabel("distance $r$")
ax.set_ylabel(r"$\langle S_z(0)S_z(r)\rangle$")
ax.set_title("TFIM (g={}, D={}): VUMPS vs iDMRG static correlator".format(g, D))
ax.legend()
fig.tight_layout()
fig.savefig("vumps_correlator_TFIM.png", dpi=150)
print("\nSaved plot to vumps_correlator_TFIM.png")
