# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator <N_0(z) N_0(0)> of an interacting fermionic chain
# with an imaginary onsite potential (gain and loss), so H is not
# Hermitian. We use the non-Hermitian KPM and compare the python
# backend with ITensor v3, which should give the same curve.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
def get_correlator(version):
    fc = fermionchain.Fermionic_Chain(n)
    if version=="python": fc.setup_python()
    else: fc.setup_cpp(version=version)
    h = 0
    for i in range(n-1):
        t = 1.0 + 0.15*i # site dependent hopping
        h = h + t*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
        h = h + 0.3*(fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    for i in range(n): # gain and loss
        h = h + 1j*0.5*(-1)**i*(1.0+0.1*i)*fc.Cdag[i]*fc.C[i]
    h = h + 0.05*fc.N[0] # lifts a degeneracy of the ground state
    fc.set_hamiltonian(h)
    es = np.linspace(0., 4., 10)
    return fc.get_dynamical_correlator(mode="DMRG", submode="KPM",
            name=[fc.N[0], fc.N[0]], E_max=12, n=60, delta=0.3, es=es)

(e1, y1) = get_correlator(3)
(e2, y2) = get_correlator("python")
print("Max difference v3-python =", np.max(np.abs(y1-y2)))

import matplotlib.pyplot as plt
plt.plot(e1, y1.real, linewidth=4, label="v3")
plt.plot(e2, y2.real, "--", linewidth=2, label="python")
plt.xlabel("energy w")
plt.ylabel("Re <N_0(z) N_0(0)>")
plt.legend()
plt.show()
