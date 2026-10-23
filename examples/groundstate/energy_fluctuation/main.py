# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Energy fluctuation sqrt(<H^2> - <H>^2) of the DMRG ground state of the
# S=1 Heisenberg chain, as a function of the bond dimension.
# It is zero only for an exact eigenstate, so it tells how
# far a variational state is from the true ground state.
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.nsweeps = 10 # number of DMRG sweeps

maxms = [1, 2, 5, 10, 20, 30, 40] # bond dimensions
des = []
for maxm in maxms:
    sc.maxm = maxm # bond dimension
    e = sc.gs_energy()
    des.append(sc.gs_energy_fluctuation())
    print("maxm", maxm, "energy", e, "fluctuation", des[-1])

import matplotlib.pyplot as plt
plt.semilogy(maxms, des, "o-")
plt.xlabel("bond dimension")
plt.ylabel("energy fluctuation")
plt.show()
