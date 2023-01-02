# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 12
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

# add non Hermitian part
for i in range(n): h = h + 1j*0.4*sc.Sz[i]

sc.set_hamiltonian(h) # initialize the system
e0 = sc.gs_energy(mode="ED") # get the ground state energy

# now compute the GS energy with MPS using different numbers of iterations 
ites = [1,4,10,40] # the default is 40, make it bigger for better accuracy

for ite in ites: # loop over number of iterations
    sc.set_hamiltonian(h) # initialize the system
    # the parameter maxit controls the number of iterations of the
    # MPS arnoldi algorithm, the bigger the more precise
    e = sc.gs_energy(mode="DMRG",maxit=ite) # get the ground state energy
    de = np.abs(e0-e) # difference with respect to the true GS energy
    print("# iterations = ",ite,", error = ",de)






