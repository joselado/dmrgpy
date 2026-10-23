# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of the S=1/2 Heisenberg chain with the Kernel
# Polynomial Method (KPM), S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0).
# We compute it with the ITensor v3 backend and with the pure Python
# backend, which must give the same spectrum.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites

def get_spectrum(version):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    h = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    name = (sc.Sz[0], sc.Sz[0])
    return sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name)

(x1, y1) = get_spectrum(3)
(x2, y2) = get_spectrum("python")
print("Max difference v3-python =", np.max(np.abs(np.array(y1)-np.array(y2))))

import matplotlib.pyplot as plt
plt.plot(x1, np.real(y1), linewidth=4, label="ITensor v3")
plt.plot(x2, np.real(y2), "--", linewidth=2, label="python")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
