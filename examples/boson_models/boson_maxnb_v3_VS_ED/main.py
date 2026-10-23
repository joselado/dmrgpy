# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of a Bose-Hubbard chain with a different number of
# boson levels maxnb on each site, H = sum t (a_i^+ a_j + h.c.) + U n (n-1).
# The backend must respect maxnb on every site, so DMRG with ITensor v3
# and with the pure Python backend must agree with exact diagonalization (ED).
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import bosonchain

n = 3
maxnb = [3, 5, 3] # levels of each site

def get_energy(U, version, mode="DMRG"):
    bc = bosonchain.Bosonic_Chain(n, maxnb=maxnb)
    if version == "python": bc.setup_python()
    else: bc.setup_cpp(version)
    bc.maxm = 120 # bond dimension
    bc.nsweeps = 50 # number of DMRG sweeps
    np.random.seed(11) # same hoppings in all runs
    h = 0
    for i in range(n-1):
        h = h + np.random.random()*(bc.Adag[i]*bc.A[i+1] + bc.Adag[i+1]*bc.A[i])
    for i in range(n):
        h = h + U*bc.N[i]*(bc.N[i]-1.0)
    bc.set_hamiltonian(h)
    return bc.gs_energy(mode=mode)

Us = [0.3, 0.7, 1.1] # onsite interaction
e3 = [get_energy(U, 3) for U in Us]
epy = [get_energy(U, "python") for U in Us]
eed = [get_energy(U, 3, mode="ED") for U in Us]
print("Max difference v3-ED =", np.max(np.abs(np.array(e3)-eed)))
print("Max difference python-ED =", np.max(np.abs(np.array(epy)-eed)))

plt.plot(Us, eed, linewidth=4, label="ED")
plt.plot(Us, e3, "--", linewidth=2, label="DMRG v3")
plt.plot(Us, epy, ":", linewidth=2, label="DMRG python")
plt.xlabel("onsite interaction U")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
