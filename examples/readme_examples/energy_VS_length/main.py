# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')


## Ground state energy of an S=1/2 spin chain, swept over the chain length
# This is the README's first snippet ("Ground state energy of an S=1/2
# spin chain", a single 30-site number) turned into the sweep the
# directory name promises: the same Heisenberg Hamiltonian is rebuilt at
# several lengths and its ground state energy computed for each. The last
# point, n=30, is exactly the README's chain, so the printed energy there
# reproduces the snippet's number.
from dmrgpy import spinchain
import numpy as np

ns = [4,8,12,16,20,24,30] # chain lengths to sweep over
es = [] # storage for the ground state energies
for n in ns: # loop over the length of the chain
  spins = ["S=1/2" for i in range(n)] # spins in each site
  sc = spinchain.Spin_Chain(spins) # create spin chain object
  h = 0 # initialize Hamiltonian
  for i in range(len(spins)-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
  sc.set_hamiltonian(h) # create the Hamiltonian
  e = sc.gs_energy() # get the ground state energy
  es.append(e) # store it
  print("Ground state energy",e,"for n =",n)

ns = np.array(ns) ; es = np.array(es)


import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 14})

# The energy itself is dominated by the trivial linear growth with n, so
# plot the energy *density* E/n: it is the quantity that converges, and
# with open boundaries it approaches the Bethe-ansatz value of the
# infinite chain, 1/4-ln(2), from above as the O(1) edge contribution is
# diluted.
plt.plot(ns,es/ns,marker="o",c="red",label="DMRG, open chain")
plt.axhline(0.25-np.log(2.),c="black",ls="--",label="Bethe ansatz, $n=\\infty$")
plt.ylabel("Energy per site $E_0/n$")
plt.xlabel("Chain length $n$")
plt.legend()

plt.tight_layout()
plt.show()
