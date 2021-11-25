# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

from dmrgpy import spinchain
n = 18 # take n sites
spins = ["S=1" for i in range(n)] # spin 1 Heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain 
sc.maxm = 20
sc.nsweeps = 30
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
# compute n excited states
es = sc.get_excited(n=6) # compute excited states with DMRG
# for the Haldane chain with open boundary conditions,
# in the limit of infnite length we expect a four-fold degeneracy
# (that correspond to the dangling S=1/2 on the edges)
# so that the bulk gap is computed as the energy difference
# between 1st and the 5th state
gap = es[4] - es[0]
print("Haldane gap =",gap)
print("Energies",es-es[0])











