from __future__ import print_function
import sys
import os
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 40 # take n sites
spins = [3 for i in range(n)] # spin 1 Heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain (open boundary)
# compute n excited states
es = sc.get_excited(n=5,mode="DMRG") # compute excited states with DMRG

# for the Haldane chain with open boundary conditions,
# in the limit of infnite length we expect a four-fold degeneracy
# (that correspond to the dangling S=1/2 on the edges)
# so that the bulk gap is computed as the energy difference
# between 1st and the 5th state

print("Excited states")
print(es)
print("\n")

sc.clean()

print("Haldane gap =",es[4]-es[0])
