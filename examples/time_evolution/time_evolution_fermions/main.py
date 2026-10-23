# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Time dependent correlator <A(t) B(0)> = <GS| A(t) B |GS> for the spinless
# fermion chain with nearest neighbor repulsion V, with A = c^dag_0 and
# B = n_0. We compute it with DMRG and compare with exact diagonalization.
import numpy as np
from dmrgpy import fermionchain, timedependent

n = 8 # number of sites
V = 4.0 # nearest neighbor repulsion
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + V*fc.N[i]*fc.N[i+1]
h = h + h.get_dagger() # hopping and interaction
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

nt, dt = 300, 0.01 # number of time steps and time step
A, B = fc.Cdag[0], fc.N[0]
(t1, y1) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="DMRG", A=A, B=B)
(t2, y2) = timedependent.evolution_ABA(fc, nt=nt, dt=dt, mode="ED", A=A, B=B)
print("Max difference DMRG-ED =", np.max(np.abs(np.array(y1)-np.array(y2))))

import matplotlib.pyplot as plt
plt.plot(t2, np.real(y2), linewidth=4, label="ED, real part")
plt.plot(t1, np.real(y1), "--", linewidth=2, label="DMRG, real part")
plt.plot(t2, np.imag(y2), linewidth=4, label="ED, imaginary part")
plt.plot(t1, np.imag(y1), "--", linewidth=2, label="DMRG, imaginary part")
plt.xlabel("time")
plt.ylabel("<A(t) B(0)>")
plt.legend()
plt.show()
