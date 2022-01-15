# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0 # initialize Hamiltonian
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

hm = sc.get_full_matrix(h) # get the operator as an sparse matrix
print("Operator as sparse matrix")
print(hm)
print("Operator as dense numpy matrix")
print(hm.todense())










