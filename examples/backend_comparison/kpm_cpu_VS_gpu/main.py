# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# KPM spectrum S(w) of Sz_0 in the S=1/2 Heisenberg chain, with the tensors
# of the pure-Python backend held on the CPU (NumPy) or on the GPU (JAX).
# Both must give the same spectrum, and the GPU only pays off for large
# KPM bond dimension. Here kpmmaxm is small, so it runs on any machine.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.pyitensor import backend

n = 8 # number of sites
es = np.linspace(-8.0, 8.0, 1000) # energies

def spectrum(name):
    backend.set_backend(name) # "numpy" or "jax"
    np.random.seed(1234) # same random start
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 30 # bond dimension
    sc.kpmmaxm = 30 # bond dimension of the KPM vectors
    return sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                    name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.15)

x1, y1 = spectrum("numpy")
x2, y2 = spectrum("jax")
print("Max difference =", np.max(np.abs(np.array(y1)-np.array(y2))))

import matplotlib.pyplot as plt
plt.plot(x1, np.real(y1), linewidth=4, label="CPU (NumPy)")
plt.plot(x2, np.real(y2), "--", linewidth=2, label="GPU (JAX)")
plt.xlim(-1, 6)
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
