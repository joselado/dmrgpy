# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # generate a Heisenberg Hamiltonian
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

h = -h

mode = "DMRG"
sc.set_hamiltonian(h)
e0 = sc.gs_energy(mode=mode) # get ground state energy
print("Energy",e0)
from dmrgpy.degeneracy import pole_eigenvalue_degeneracy

# sweep the pole broadening "delta" used by the degeneracy estimator;
# this is the cheap parameter available in this example and the
# degeneracy estimate should stay close to the true value (3, the
# ferromagnetic Heisenberg dimer's triplet ground state) across it
deltas = [0.05,0.1,0.2,0.3,0.5]
degs = []
for delta in deltas:
    deg = pole_eigenvalue_degeneracy(sc,h,e0,mode=mode,delta=delta) # compute the degeneracy
    print("Delta",delta,"Degeneracy",deg)
    degs.append(deg.real)

print("Excited states",sc.get_excited(n=int(degs[1])+2,mode=mode))

import matplotlib.pyplot as plt
plt.plot(deltas,degs,marker="o",label="estimated degeneracy")
plt.axhline(3.0,color="gray",linestyle="--",label="expected (3)")
plt.xlabel("Pole broadening delta")
plt.ylabel("Ground state degeneracy")
plt.legend()
plt.show()
