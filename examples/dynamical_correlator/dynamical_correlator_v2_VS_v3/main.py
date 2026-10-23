# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0)
# of a Heisenberg chain of S=1 spins, computed with the KPM-DMRG method.
# We run the two ITensor backends (v2 and v3) and compare them.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
es = np.linspace(-0.5, 4.5, 200) # energies
results = []
for version in [2, 3]:
    sc = spinchain.Spin_Chain(["S=1" for i in range(n)], itensor_version=version)
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
    results.append(sc.get_dynamical_correlator(name=name, es=es, delta=0.1,
        mode="DMRG"))

(e2, s2), (e3, s3) = results
print("Max difference v2-v3 =", np.max(np.abs(s2 - s3)))

import matplotlib.pyplot as plt
plt.plot(e2, s2.real, linewidth=4, label="ITensor v2")
plt.plot(e3, s3.real, "--", linewidth=2, label="ITensor v3")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
