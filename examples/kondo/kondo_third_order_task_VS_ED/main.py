# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Third-order STM/Kondo spectrum dI/dV at T=0 of an S=1 impurity with
# easy-plane anisotropy in a field, coupled to a chain of S=1/2 spins.
# We use DMRG with submode="TASK" (tangent-space Krylov) and compare it
# with exact diagonalization. A small maxm moves the curve.
import numpy as np
from dmrgpy import spinchain

G, MUB = 2.0, 5.7883818066e-5 # g factor, Bohr magneton in eV/T
nb = 6 # number of S=1/2 spins

def chain(maxm):
    sc = spinchain.Spin_Chain(["S=1"]+["S=1/2"]*nb, itensor_version="python")
    h = 1e-3*sc.Sz[0]*sc.Sz[0] + G*MUB*2.0*sc.Sz[0] # D=1 meV, B=2 T
    for i in range(nb):
        J = 0.6e-3 if i==0 else 1e-3 # exchange coupling to the impurity
        h = h + J*(sc.Sx[i]*sc.Sx[i+1]+sc.Sy[i]*sc.Sy[i+1]+sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    sc.maxm = maxm # bond dimension
    sc.nsweeps = 20 # number of DMRG sweeps
    sc.cutoff = 0.
    return sc

eV = np.linspace(-5e-3, 5e-3, 201) # bias voltage
kw = dict(site=0, Jrho_s=-0.1, U=0.3, T=0., omega0=20e-3, Gamma0=5e-5, order=3)
sc = chain(64) # large enough to hold the whole Hilbert space
_, ed = sc.get_kondo_spectrum(eV, mode="ED", **kw)
_, t1 = sc.get_kondo_spectrum(eV, mode="DMRG", submode="TASK", nkry=200, **kw)
_, t2 = chain(4).get_kondo_spectrum(eV, mode="DMRG", submode="TASK", nkry=200, **kw)
print("Max difference =", np.max(np.abs(t1-ed))/np.max(np.abs(ed)))
assert np.max(np.abs(t1-ed)) < 1e-8*np.max(np.abs(ed))

import matplotlib.pyplot as plt
plt.plot(eV*1e3, ed, "k", linewidth=4, label="ED")
plt.plot(eV*1e3, t1, "--", linewidth=2, label="TASK, maxm=64")
plt.plot(eV*1e3, t2, ":", linewidth=2, label="TASK, maxm=4")
plt.xlabel("eV (meV)")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
