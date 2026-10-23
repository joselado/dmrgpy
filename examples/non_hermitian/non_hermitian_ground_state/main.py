# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Selecting a total spin sector with a non-Hermitian term. We add the
# imaginary term i(S^2 - S0(S0+1)) to the Heisenberg chain, with S^2 the
# total spin operator. Its eigenvalues stay real only for S=S0. We ask
# get_gs() for the state it picks and measure its total spin S(S+1).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 10 # bond dimension

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
s2 = sum(sc.Sx)*sum(sc.Sx) + sum(sc.Sy)*sum(sc.Sy) + sum(sc.Sz)*sum(sc.Sz)

targets = [0., 1., 2.] # total spin S0
found = []
for S0 in targets:
    sc.set_hamiltonian(h + 1j*(s2 - S0*(S0+1.))) # non-Hermitian Hamiltonian
    wf = sc.get_gs()
    found.append(wf.dot(s2*wf).real) # S(S+1) of the state
    print("Target S0 =", S0, " S(S+1) =", found[-1])

import matplotlib.pyplot as plt
x = [S0*(S0+1.) for S0 in targets]
plt.plot(x, x, "k--", label="target")
plt.plot(x, found, "o", markersize=10, label="get_gs()")
plt.xlabel("target S0(S0+1)")
plt.ylabel("S(S+1) of the state")
plt.legend()
plt.show()
