# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Thermal expectation value <Sz_0> = Tr(exp(-H/T) Sz_0)/Z of a Heisenberg
# chain in a magnetic field B. METTS samples minimally entangled thermal
# states, so it works with MPS at finite temperature. It is a Monte Carlo
# method, so the error bars are statistical. We check it against ED.
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
B = 0.3 # magnetic field
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")

# Heisenberg Hamiltonian plus a field, H = sum_i S_i . S_{i+1} + B sum_i Sz_i
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + B*sc.Sz[i]
sc.set_hamiltonian(h)

Ts = [0.5, 1.0, 2.0, 3.0] # temperatures
ed, metts, err = [], [], []
for T in Ts:
    ed.append(sc.vev(sc.Sz[0], mode="ED", T=T).real)
    (m, e) = sc.metts_vev(sc.Sz[0], T, nsamples=300, nwarmup=30, seed=42)
    metts.append(m.real) ; err.append(e)
print("Max difference METTS-ED =", np.max(np.abs(np.array(metts)-ed)))

import matplotlib.pyplot as plt
plt.plot(Ts, ed, linewidth=4, label="ED")
plt.errorbar(Ts, metts, yerr=err, fmt="o", capsize=3, label="METTS")
plt.xlabel("temperature T")
plt.ylabel("<Sz_0>")
plt.legend()
plt.show()
