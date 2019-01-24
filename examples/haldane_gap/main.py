# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

from dmrgpy import spinchain
n = 60 # take n sites
spins = [3 for i in range(n)] # spin 1 Heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain 
def fj(i,j): # function for the coupling (open boundary)
  if abs(j-i)==1: return 1.0 # first neighbor to the right
  else: return 0.0 # anything else
sc.set_exchange(fj) # Add coupling between the spins
# compute n excited states
es = sc.get_excited(n=5) # compute excited states with DMRG
# for the Haldane chain with open boundary conditions,
# in the limit of infnite length we expect a four-fold degeneracy
# (that correspond to the dangling S=1/2 on the edges)
# so that the bulk gap is computed as the energy difference
# between 1st and the 5th state
gap = es[4] - es[0]
print("Haldane gap =",gap)
sc.clean() # clean files


