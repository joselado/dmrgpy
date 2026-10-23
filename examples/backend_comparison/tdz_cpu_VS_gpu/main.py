# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator S(w) of the first site of a Heisenberg chain with
# complex time evolution (submode="TDZ"), using the pure Python backend.
# It is computed with NumPy on the CPU and with JAX, which runs on the
# GPU when there is one. Both must give the same spectrum.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.pyitensor import backend

def spectrum(name):
    backend.set_backend(name) # "numpy" or "jax"
    n = 8 # number of sites
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")
    # Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 16 # bond dimension
    sc.nsweeps = 3 # number of DMRG sweeps
    es = np.linspace(-0.5, 4.0, 200) # energies
    return sc.get_dynamical_correlator(mode="DMRG", submode="TDZ", es=es,
            name=(sc.Sz[0], sc.Sz[0]), dt=0.1, nt=20)

np.random.seed(1234) # same random start for both runs
(e1, s1) = spectrum("numpy")
np.random.seed(1234)
(e2, s2) = spectrum("jax")
backend.set_backend("numpy")
print("Max difference CPU-GPU =", np.max(np.abs(np.array(s1)-np.array(s2))))

import matplotlib.pyplot as plt
plt.plot(e1, np.real(s1), linewidth=4, label="CPU (NumPy)")
plt.plot(e2, np.real(s2), "--", linewidth=2, label="GPU (JAX)")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
