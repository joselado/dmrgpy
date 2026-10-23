# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Time evolution of a state after applying an operator A to it,
# <X|A^dag e^{iHt} B e^{-iHt} A|X> = <psi(t)|B|psi(t)>, with
# |psi(t)> = e^{-iHt} A|X>. Here |X> is the Neel state of a 3 site
# Heisenberg chain, A = Sx_0 + Sy_1 and B = Sz_0. We use DMRG and ED.
import numpy as np
from dmrgpy import spinchain, timedependent

n = 3 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h0 = 0 # staggered field, its ground state is the Neel state
for i in range(n): h0 = h0 + (-1)**i*sc.Sz[i]
h1 = 0 # Heisenberg Hamiltonian
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
wf, wfED = sc.get_gs(), sc.get_gs(mode="ED") # Neel states
sc.set_hamiltonian(h1) # Hamiltonian for the evolution

A = sc.Sx[0] + sc.Sy[1] # operator applied to the state
B = sc.Sz[0] # operator that we measure
nt, dt = 500, 0.02 # number of time steps, time step
(ts, sz) = timedependent.evolution_ABA(sc, A=A, B=B, nt=nt, dt=dt, wf=wf)
(tsED, szED) = timedependent.evolution_ABA(sc, A=A, B=B, nt=nt, dt=dt,
                                           wf=wfED, mode="ED")
print("Max difference =", np.max(np.abs(sz-szED)))

import matplotlib.pyplot as plt
plt.plot(tsED, szED.real, linewidth=4, label="ED")
plt.plot(ts, sz.real, "--", linewidth=2, label="DMRG")
plt.xlabel("time t")
plt.ylabel("<B>(t)")
plt.legend()
plt.show()
