# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 3 # number of sites in the spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
Sx,Sy,Sz = sc.Sx,sc.Sy,sc.Sz # lists with the respective operators

# generate a Heisenberg operator as example
A0 = 0
for i in range(n-1):
    A0 = A0 + Sx[i]*Sx[i+1]
    A0 = A0 + Sy[i]*Sy[i+1]
    A0 = A0 + Sz[i]*Sz[i+1]

for i in range(n): # add some terms that make the operator non-Hermitian
    A0 = A0 + 0.3*1j*Sy[i]*Sx[i]
A0 = A0 + 0.2*1j # add a diagonal imaginary
# this controls the accuracy of the tensor-network calculation
# sc.mpomaxm = 80 
# sc.maxm = 80 
# (the higher, the more accurate) typical minimal values are 50-100


# in order to make products in a scalable way, convert to a new object
A = sc.toMPO(A0) # transforms into an MPO (that can be multiplied)
# the object A is a "StaticOperator" you can do A*B  but not A + B
#A = sc.toMPO(A0,mode="ED") # uncomment to test it with ED using this instead

Ad = A.get_dagger() # compute the Hermitian

# define a few operators via multiplying
AA = A*A
AdA = Ad*A
AAd = A*Ad
AAdA = A*Ad*A

print("Trace of A",A.trace())
print("Trace of A^\dagger",Ad.trace())
print("Trace of A*A",AA.trace())
print("Trace of A*A^\dagger",AAd.trace())
print("Trace of A*A^\dagger*A",AAdA.trace())


