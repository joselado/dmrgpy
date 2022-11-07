# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6 # number of sites in the spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
Sx,Sy,Sz = sc.Sx,sc.Sy,sc.Sz # lists with the respective operators

# generate a Heisenberg model Hamiltonian as example
h = 0
for i in range(n-1):
    h = h + Sx[i]*Sx[i+1]
    h = h + Sy[i]*Sy[i+1]
    h = h + Sz[i]*Sz[i+1]

# this controls the accuracy of the tensor-network calculation
# (the higher, the more accurate) typical minimal values are 50-100
sc.maxm = 200 

Op = sc.toMPO(h) # transforms into an MPO (that can be multiplied)
#Op = sc.toMPO(h,mode="ED") # uncomment to test it with ED using this instead

## Op objects have two operations, product and trace
# Op_C = Op_A*Op_B this is equivalent to the matrix-product of operators
# c = Op.trace() this is the trace of the operator


N = 4 # compute some products
Opk = Op.copy() # make a copy
for i in range(N):
    out = Sz[0]*Opk*Sz[1] # compute the product
    Opk = Opk*Op # compute the next power
    print("Trace of S^z_0 H^"+str(i)+" S^z_N",out.trace())






