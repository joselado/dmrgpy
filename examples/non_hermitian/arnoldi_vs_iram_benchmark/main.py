# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of a non-Hermitian interacting fermion chain, with gain and
# loss i(-1)^i on the sites. The ground state is the eigenvalue with the
# smallest real part. We compare two Arnoldi solvers: the thick restart one
# (arnolditk) and the implicitly restarted one (arpacktk), by the number of
# H*x products they need and by the error of the energy against ED.
import numpy as np
from dmrgpy import fermionchain
from dmrgpy.edtk.edchain import EDOperator, State
from dmrgpy.algebra import arnolditk, arpacktk

count = [0] # number of products H*x
mul = EDOperator.__mul__
def counted_mul(self, v):
    if isinstance(v, State): count[0] += 1
    return mul(self, v)
EDOperator.__mul__ = counted_mul
arnolditk.arnoldimode = "ED"

ns = [4, 6, 8] # number of sites
nops, err = {"arnoldi": [], "iram": []}, {"arnoldi": [], "iram": []}
for n in ns:
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    for i in range(n): # gain and loss
        h = h + 1j*(-1)**i*fc.N[i]
    for i in range(n-1): # interaction
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    es = np.array(fc.get_excited(mode="ED", n=6)) # exact energies
    for name in nops:
        count[0] = 0
        if name == "arnoldi": e = arnolditk.mpsarnoldi(fc, h, mode="GS", delta=1e-6)[0]
        else: e = arpacktk.mpsiram(fc, h, which="SR", tol=1e-6)[0]
        nops[name].append(count[0])
        err[name].append(max(np.min(np.abs(np.atleast_1d(e)[0]-es)), 1e-16))
print("H*x products:", nops)
print("Energy errors:", err)

import matplotlib.pyplot as plt
for name in nops:
    plt.plot(ns, nops[name], marker="o", label=name)
plt.xlabel("number of sites")
plt.ylabel("number of H*x products")
plt.legend()
plt.show()
