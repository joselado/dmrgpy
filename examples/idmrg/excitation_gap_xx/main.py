# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
import matplotlib.pyplot as plt  # plotting
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

###########################################################################
### First excited state of an infinite chain: the tangent-space/        ###
### quasiparticle excitation ansatz (pyitensor/idmrg_excitations.py)    ###
###########################################################################
# Unlike a finite chain, an infinite chain's excitations form a momentum-
# resolved band E(k), not a single discrete state -- see
# idmrg_excitations.py's own module docstring for the algorithm. This
# feature requires gs_method="vumps" (not the default "idmrg") -- see
# infinitechain.py's own class docstring for why.
#
# This example uses a field-polarized XX chain (field above the model's
# own saturation field), whose ground state is exactly the fully-polarized
# (D=1) product state, letting the whole dispersion be checked against the
# exact (free-fermion) answer; see excitation_gap_tfim for a genuinely
# entangled (D>1) example.

J, h = 1.0, 3.0  # h > J: above the XX chain's own saturation field
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "vumps"
H = -J*(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]) - 2*h*ic.SzC[0]
ic.set_hamiltonian(H)

ic.maxm = 4     # D=1 ground state here, this is just a safety cap
ic.maxiter = 300
ic.etol = 1e-12

e0 = ic.gs_energy()
print("VUMPS converged:", ic.converged)
print("ground-state energy density:", e0, "(exact: -{})".format(h))
print()

ks = np.linspace(0, 2*np.pi, 13)
exact = 2*h - J*np.cos(ks)
got = np.array([ic.excitation_energies(k, n=1)[0].real for k in ks])
print("momentum-resolved dispersion E(k) (exact: 2h - J*cos(k)):")
for k, e_k, e_exact in zip(ks, got, exact):
    print("  k={:.4f}  E(k)={:.6f}  exact={:.6f}  diff={:.2e}".format(
        k, e_k, e_exact, e_k - e_exact))
assert np.allclose(got, exact, atol=1e-6)

gap = ic.excitation_gap()
print()
print("gap = min_k E(k):", gap, " (exact: 2h - J =", 2*h - J, ")")
assert abs(gap - (2*h - J)) < 1e-6

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(ks, exact, "-", color="black", label="exact (free fermion)")
ax.plot(ks, got, "o", color="tab:red", markersize=5, label="dmrgpy (VUMPS D=1)")
ax.set_xlabel("momentum $k$")
ax.set_ylabel("$E(k)$")
ax.set_title("Field-polarized XX chain: single-magnon dispersion")
ax.legend()
fig.tight_layout()
fig.savefig("excitation_gap_xx.png", dpi=150)
print("\nsaved plot to excitation_gap_xx.png")
