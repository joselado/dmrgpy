# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain

def get_energy(n):
    spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(n-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc.gs_energy() # compute the ground state energy

n = 20
e = get_energy(n)
print("Energy",e)

# sweep the chain length and look at the ground state energy per bond,
# which converges to the known Bethe-ansatz value (1/4 - ln(2)) for the
# infinite spin-1/2 Heisenberg chain
ns = [4,8,12,16,20]
es = [get_energy(ni) for ni in ns]
eperbond = [e/(ni-1) for e,ni in zip(es,ns)]
ebethe = 0.25 - np.log(2)

plt.plot(ns,eperbond,marker="o",label="DMRG")
plt.axhline(ebethe,color="black",linestyle="--",label="Bethe ansatz (n -> inf)")
plt.xlabel("Chain length")
plt.ylabel("Ground state energy per bond")
plt.legend()
plt.show()
