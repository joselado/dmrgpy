# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench of a spin chain. We prepare the ground state of a staggered
# field, H0 = sum_i (-1)^i Sz_i, and evolve it with the Heisenberg
# Hamiltonian, <Sz_0>(t) = <GS|exp(iHt) Sz_0 exp(-iHt)|GS>.
# This is the dynamics seen in cold atom quench experiments. We compare
# the pure Python TDVP with ITensor v3.
import numpy as np
from dmrgpy import spinchain, timedependent

n = 6 # number of sites

def evolve(version):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version=version)
    h0 = 0
    for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i] # staggered field
    h1 = 0
    for i in range(n-1): # Heisenberg Hamiltonian
        h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.maxm = 40 # bond dimension
    sc.nsweeps = 10 # number of DMRG sweeps
    sc.set_hamiltonian(h0)
    wf = sc.get_gs() # initial state
    sc.set_hamiltonian(h1) # quench
    return timedependent.evolve_and_measure(sc, operator=sc.Sz[0], nt=200, dt=0.02, wf=wf)

(t3, s3) = evolve(3)
(tp, sp) = evolve("python")
print("Max difference =", np.max(np.abs(np.array(s3)-np.array(sp))))

import matplotlib.pyplot as plt
plt.plot(t3, np.real(s3), linewidth=4, label="v3")
plt.plot(tp, np.real(sp), "--", linewidth=2, label="python")
plt.xlabel("time t")
plt.ylabel("<Sz_0>")
plt.legend()
plt.show()
