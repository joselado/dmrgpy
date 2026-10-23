# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Gap of the transverse-field Ising chain, H = sum_i Sz_i Sz_{i+1} + B sum_i Sx_i.
# The gap is the difference between the two lowest energies, and it is
# what a spectroscopy measurement sees first. We compute the excited
# states with ITensor v3 and with the pure Python backend, as a function of B.
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites

def get_gap(B, version):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    h = 0
    for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1]
    for i in range(n): h = h + B*sc.Sx[i]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    es = sc.get_excited(mode="DMRG", n=2)
    return es[1] - es[0]

Bs = np.linspace(0.2, 1.0, 9) # transverse fields
g3 = [get_gap(B, 3) for B in Bs]
gp = [get_gap(B, "python") for B in Bs]
print("Max difference v3-python =", np.max(np.abs(np.array(g3)-gp)))

import matplotlib.pyplot as plt
plt.plot(Bs, g3, linewidth=4, label="ITensor v3")
plt.plot(Bs, gp, "--", linewidth=2, label="python")
plt.xlabel("transverse field B")
plt.ylabel("gap")
plt.legend()
plt.show()
