# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlator <Sz_0 Sz_r> in the ground state of the infinite transverse
# field Ising chain, H = -sum_i sx_i sx_{i+1} - g sum_i sz_i, in the gapped
# phase g>1. The correlation decays exponentially with the distance r.
# We use the VUMPS ground state and compare with the iDMRG ground state.
import numpy as np
from dmrgpy import infinitechain

g = 1.5 # transverse field

def get_chain(gs_method):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
    ic.set_hamiltonian(-4*ic.SxC[0]*ic.SxR[0] - 2*g*ic.SzC[0]) # Ising + field
    ic.gs_method = gs_method
    ic.maxm = 3 # bond dimension
    return ic

ic1 = get_chain("vumps")
ic1.maxiter = 800
ic1.vumps_nrestarts = 6
ic2 = get_chain("idmrg")
ic2.maxiter = 200
ic2.etol = 1e-10
print("Energy per site: VUMPS =", ic1.gs_energy(), "iDMRG =", ic2.gs_energy())

rs = range(8) # distances
c1 = np.array([ic1.correlator("Sz", 0, "Sz", r).real for r in rs])
c2 = np.array([ic2.correlator("Sz", 0, "Sz", r).real for r in rs])
print("Max difference VUMPS-iDMRG =", np.max(np.abs(c1-c2)))
assert np.max(np.abs(c1-c2)) < 0.1

import matplotlib.pyplot as plt
plt.plot(rs, c1, linewidth=4, label="VUMPS")
plt.plot(rs, c2, "--", linewidth=2, label="iDMRG")
plt.xlabel("distance r")
plt.ylabel("<Sz_0 Sz_r>")
plt.legend()
plt.show()
