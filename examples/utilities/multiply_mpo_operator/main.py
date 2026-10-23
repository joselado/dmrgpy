# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Operators can be written as matrix product operators (MPO) and
# multiplied. Here we build the Heisenberg Hamiltonian H as an MPO and
# compute the trace of Sz_0 H^k Sz_1 for several powers k.
# Traces like this are the moments of a spectral function.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.maxm = 200 # bond dimension of the MPO

Op = sc.toMPO(h) # Hamiltonian as an MPO
Opk = Op.copy() # this will be the k-th power of H
traces = []
for k in range(1, 5):
    out = sc.Sz[0]*Opk*sc.Sz[1] # product of three operators
    traces.append(out.trace())
    Opk = Opk*Op # next power of H
    print("k =", k, "trace =", traces[-1])

import matplotlib.pyplot as plt
plt.plot(range(1, 5), np.real(traces), marker="o", label="real part")
plt.plot(range(1, 5), np.imag(traces), marker="o", label="imaginary part")
plt.xlabel("power k of H")
plt.ylabel("trace of Sz_0 H^k Sz_1")
plt.legend()
plt.show()
