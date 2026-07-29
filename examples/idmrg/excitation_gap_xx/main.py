# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np  # conventional numpy library
from dmrgpy import infinitechain  # infinite-DMRG (iDMRG) chain object

###########################################################################
### First excited state of an infinite chain: the tangent-space/        ###
### quasiparticle excitation ansatz (pyitensor/idmrg_excitations.py)    ###
###########################################################################
# Unlike a finite chain, an infinite chain's excitations form a momentum-
# resolved band E(k), not a single discrete state -- see
# idmrg_excitations.py's own module docstring for the algorithm.
#
# IMPORTANT SCOPE NOTE: only a product-state-like (bond dimension D=1)
# converged ground state is currently supported -- see idmrg_excitations.py's
# "KNOWN LIMITATION" section for what was tried (and ruled out) trying to
# extend this to a genuinely entangled (D>1) ground state. This example
# therefore uses a field-polarized XX chain (field above the model's own
# saturation field), whose ground state is exactly the fully-polarized
# product state, so iDMRG converges to D=1 -- matching the currently
# supported scope, and letting the whole dispersion be checked against the
# exact (free-fermion) answer.

J, h = 1.0, 3.0  # h > J: above the XX chain's own saturation field
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
H = -J*(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]) - 2*h*ic.SzC[0]
ic.set_hamiltonian(H)

ic.maxm = 4     # D=1 ground state here, this is just a safety cap
ic.maxiter = 60
ic.etol = 1e-12

e0 = ic.gs_energy()
print("iDMRG converged:", ic.converged)
print("ground-state energy density:", e0, "(exact: -{})".format(h))
print()

ks = np.linspace(0, 2*np.pi, 13)
print("momentum-resolved dispersion E(k) (exact: 2h - J*cos(k)):")
for k in ks:
    e_k = ic.excitation_energies(k, n=1)[0]
    exact = 2*h - J*np.cos(k)
    print("  k={:.4f}  E(k)={:.6f}  exact={:.6f}  diff={:.2e}".format(
        k, e_k, exact, e_k - exact))

print()
gap = ic.excitation_gap()
print("gap = min_k E(k):", gap, " (exact: 2h - J =", 2*h - J, ")")
