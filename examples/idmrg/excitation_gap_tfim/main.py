# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt  # plotting
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

###########################################################################
### First excited state of an infinite chain, D>1 (genuinely entangled) ###
### ground state: the tangent-space/quasiparticle excitation ansatz     ###
### (pyitensor/idmrg_excitations.py), transverse-field Ising model      ###
###########################################################################
# Unlike excitation_gap_xx (a D=1, exactly-product-state example), the
# transverse-field Ising ground state below is genuinely entangled at
# D=2 -- see idmrg_excitations.py's own module docstring/"History"
# section for the tangent-space excitation ansatz that makes this D>1
# case exact (matching MPSKit.jl's own independent D=2 result to 6
# significant figures, see docs/idmrg_excitation_mpskit_port_plan.md).
# Requires gs_method="vumps" -- see infinitechain.py's own class
# docstring for why.
#
# H = -4*SxC[i]*SxC[i+1] - 2*g*SzC[i] = -sigma^x_i sigma^x_{i+1} - g*sigma^z_i
# (Sx=sigma^x/2, Sz=sigma^z/2), the standard Pauli-matrix TFIM convention.
# Exact free-fermion single-magnon dispersion:
# eps(k) = 2*sqrt(J^2 + g^2 - 2*J*g*cos(k)), J=1.

J, g = 1.0, 2.5  # g > 1: gapped paramagnetic phase
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "vumps"
H = -4.0*J*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0]
ic.set_hamiltonian(H)

ic.maxm = 2   # target VUMPS bond dimension D=2
ic.etol = 1e-12
ic.vumps_nrestarts = 6

e0 = ic.gs_energy()
print("VUMPS converged:", ic.converged)
print("ground-state energy density:", e0)
print()

ks = np.linspace(0, np.pi, 13)
exact = 2*np.sqrt(J**2 + g**2 - 2*J*g*np.cos(ks))
got = np.array([ic.excitation_energies(k, n=1)[0].real for k in ks])
print("momentum-resolved dispersion E(k) (exact free-fermion band):")
for k, e_k, e_exact in zip(ks, got, exact):
    print("  k={:.4f}  E(k)={:.6f}  exact={:.6f}  diff={:.2e}".format(
        k, e_k, e_exact, e_k - e_exact))
assert np.allclose(got, exact, atol=2e-3)

gap = ic.excitation_gap(ks=np.linspace(0, np.pi, 41))
print()
print("gap = min_k E(k):", gap, " (exact:", 2*(g - J), ")")
assert abs(gap - 2*(g - J)) < 2e-3

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(ks, exact, "-", color="black", label="exact (free fermion)")
ax.plot(ks, got, "o", color="tab:red", markersize=5, label="dmrgpy (VUMPS D=2)")
ax.set_xlabel("momentum $k$")
ax.set_ylabel("$E(k)$")
ax.set_title("Transverse-field Ising (g={}): single-magnon dispersion".format(g))
ax.legend()
fig.tight_layout()
fig.savefig("excitation_gap_tfim.png", dpi=150)
print("\nsaved plot to excitation_gap_tfim.png")
