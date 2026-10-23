# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin correlator S(w) = sum_n <GS|Sz_i|n><n|Sz_j|GS> delta(w - E_n + E_0).
# We compute it with the spin KPM routine and with the general KPM
# routine, which takes any two states and the operator X that is expanded.
# Both must give the same spectrum, which is what inelastic neutrons measure.
import numpy as np
from dmrgpy import spinchain
from dmrgpy.kpmdmrg import kpm_wfa_wfb

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension in the KPM
wf = sc.get_gs() # ground state

es = np.linspace(-0.5, 6, 2000) # energies
delta = 0.05 # broadening of the peaks
(x1, y1) = sc.get_dynamical_correlator(mode="DMRG", name=(sc.Sz[0], sc.Sz[1]), es=es, delta=delta)

# general KPM, from the two states Sz_0|GS> and Sz_1|GS>
(x2, y2) = kpm_wfa_wfb(sc, wfa=sc.Sz[0]*wf, wfb=sc.Sz[1]*wf, delta=delta, X=h)
x2 = x2 - sc.gs_energy() # shift by the ground state energy

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="spin KPM")
plt.plot(x2, y2.real, "--", linewidth=2, label="general KPM")
plt.xlim([-0.5, 4.5])
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
