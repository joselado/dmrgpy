# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Non-Hermitian DMRG on the Julia backend for a Hatano-Nelson chain:
# asymmetric hopping (t_right=1, t_left) plus a staggered imaginary
# potential. The eigenvalue with the smallest real part is compared with
# exact diagonalization (ED), as a function of the hopping asymmetry.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
gamma = 0.7 # staggered imaginary potential

def chain(t_left, julia):
    fc = fermionchain.Fermionic_Chain(n)
    if julia: fc.setup_julia()
    else: fc.setup_python()
    h = 0
    for i in range(n-1): # asymmetric hopping
        h = h + fc.Cdag[i]*fc.C[i+1] + t_left*fc.Cdag[i+1]*fc.C[i]
    for i in range(n): # staggered imaginary potential
        h = h + 1j*gamma*(-1)**i*fc.Cdag[i]*fc.C[i]
    for i in range(n-1): # interaction
        h = h + 0.5*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    return fc

t_lefts = [1.0, 0.8, 0.6, 0.4, 0.2]
e_ed, e_dmrg = [], []
for t_left in t_lefts:
    e_ed.append(chain(t_left, False).get_excited(mode="ED", n=4)[0])
    e, psil, psir = chain(t_left, True).nhdmrg() # left and right states
    e_dmrg.append(e)
diff = np.max(np.abs(np.real(e_ed) - np.real(e_dmrg)))
print("Max difference =", diff)
assert diff < 1e-6

import matplotlib.pyplot as plt
plt.plot(t_lefts, np.real(e_ed), linewidth=4, label="ED")
plt.plot(t_lefts, np.real(e_dmrg), "--", linewidth=2, label="NH-DMRG")
plt.xlabel("t_left")
plt.ylabel("Re(eigenvalue)")
plt.legend()
plt.show()
