# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
from dmrgpy import spinfermionchain
n = 7
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
def geth(sc):
  h = 0
  for i in range(n-1):
      h = h + sc.Sx[i]*sc.Sx[i+1]
      h = h + sc.Sy[i]*sc.Sy[i+1]
      h = h + sc.Sz[i]*sc.Sz[i+1]
  return h

sc = spinchain.Spin_Chain(spins) # create the spin chain
fc = spinfermionchain.Spin_Fermion_Hamiltonian(["S" for s in spins])

sc.set_hamiltonian(geth(sc)) # create Hamiltonian
fc.set_hamiltonian(geth(fc)) # create Hamiltonian

print("Energy with spins",sc.gs_energy())
print("Energy with fermions",fc.gs_energy())
print(fc.get_density())

# sweep the chain length to confirm the two representations keep
# agreeing on the ground-state energy (cheap, small chains)
ns_list = [3,4,5,6,7]
e_spin = []
e_fermion = []
for nn in ns_list:
    n = nn
    spins_n = [2 for i in range(n)]
    scn = spinchain.Spin_Chain(spins_n)
    fcn = spinfermionchain.Spin_Fermion_Hamiltonian(["S" for s in spins_n])
    scn.set_hamiltonian(geth(scn))
    fcn.set_hamiltonian(geth(fcn))
    e_spin.append(scn.gs_energy())
    e_fermion.append(fcn.gs_energy())
    print("n =",nn,"E(spins) =",e_spin[-1],"E(fermions) =",e_fermion[-1])

import matplotlib.pyplot as plt
plt.plot(ns_list,e_spin,marker="o",label="Spin chain")
plt.plot(ns_list,e_fermion,marker="x",label="Spin-fermion chain")
plt.xlabel("Number of sites")
plt.ylabel("Ground-state energy")
plt.legend()
plt.show()











