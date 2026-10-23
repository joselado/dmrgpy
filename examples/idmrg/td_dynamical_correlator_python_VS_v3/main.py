# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical spin correlator S(k,w) of the infinite S=1/2 Heisenberg chain,
# from the real-time evolution of a window of the chain (iDMRG ground state).
# It is the quantity that inelastic neutron scattering measures.
# We compare the python and the v3 backends, which must give the same S(k,w).
import numpy as np
from dmrgpy import infinitechain

def get_skw(version):
    ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"], itensor_version=version)
    ic.gs_method = "idmrg"
    h = 0
    h = (ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1]
         + ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0])
    ic.set_hamiltonian(h)
    ic.maxm = 10 # bond dimension
    ic.maxiter = 40
    ic.etol = 1e-11
    ic.gs_energy()
    return ic.td_dynamical_correlator("Sz", 0, "Sz", 6, dt=0.05, nt=4,
            x_values=[-2, -1, 0, 1, 2], maxdim=20)

(ks, ws, s_py) = get_skw("python")
(ks, ws, s_v3) = get_skw(3)
print("Max difference =", np.max(np.abs(s_py-s_v3)))
ik = np.argmin(np.abs(np.abs(ks)-np.pi/2)) # momentum close to pi/2

import matplotlib.pyplot as plt
plt.plot(ws, np.abs(s_py[ik]), linewidth=4, label="python")
plt.plot(ws, np.abs(s_v3[ik]), "--", linewidth=2, label="v3")
plt.xlabel("energy w")
plt.ylabel("|S(k,w)| at k = pi/2")
plt.legend()
plt.show()
