# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlator <Sz_0 Sz_r> in the ground state of the infinite transverse-field
# Ising chain, H = -sum_i sx_i sx_{i+1} - g sum_i sz_i, in the gapped phase g>1.
# The ground state is found with VUMPS, and the correlator decays
# exponentially with the distance. We compare the pure-Python backend
# with itensor_version=3, for the bond dimension D.
import numpy as np
from dmrgpy import infinitechain

g = 1.5 # transverse field
D = 3 # bond dimension
rs = range(6) # distances

def correlator(version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=version)
    ic.gs_method = "vumps"
    ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.maxm = D
    ic.maxiter = 800 # maximum number of iterations
    ic.vumps_nrestarts = 6
    ic.gs_energy()
    return np.array([ic.correlator("Sz", 0, "Sz", r).real for r in rs])

c1 = correlator("python")
c2 = correlator(3)
print("Max difference =", np.max(np.abs(c1-c2)))
assert np.max(np.abs(c1-c2)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(rs, c2, "o", markersize=12, label="itensor_version=3")
plt.plot(rs, c1, "x", markersize=10, label="python")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.legend()
plt.show()
