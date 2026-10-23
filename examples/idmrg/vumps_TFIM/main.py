# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy density of the infinite transverse-field Ising chain,
# H = -sum sigma^x_i sigma^x_{i+1} - g sum sigma^z_i, in the gapped phase g>1.
# We compare two algorithms, VUMPS and iDMRG, as a function of the bond
# dimension D, with the exact free-fermion result.
import numpy as np
from scipy.integrate import quad
from dmrgpy import infinitechain

g = 1.5 # transverse field
exact = -quad(lambda k: np.sqrt(1 + g**2 - 2*g*np.cos(k)), 0, np.pi)[0]/np.pi

def energy(method, D):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.gs_method = method
    ic.maxm = D # bond dimension
    ic.maxiter = 400
    if method == "vumps": ic.vumps_nrestarts = 6
    else: ic.etol = 1e-10
    return ic.gs_energy()

Ds = [1, 2, 3, 4] # bond dimensions
ev = [energy("vumps", D) for D in Ds]
ei = [energy("idmrg", D) for D in Ds]
print("Max difference VUMPS-exact =", np.max(np.abs(np.array(ev)-exact)))
assert np.max(np.abs(np.array(ev)-exact)) < 0.5

import matplotlib.pyplot as plt
plt.plot(Ds, ev, "o-", label="VUMPS")
plt.plot(Ds, ei, "s-", label="iDMRG")
plt.axhline(exact, color="gray", linestyle="--", label="exact")
plt.xlabel("bond dimension D")
plt.ylabel("energy density")
plt.legend()
plt.show()
