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

mode = "ED"
sc.set_hamiltonian(h)
e0 = sc.gs_energy(mode=mode) # get ground state energy
print("Energy",e0)
from dmrgpy.degeneracy import pole_eigenvalue_degeneracy
deg = pole_eigenvalue_degeneracy(sc,h,e0,mode=mode,delta=0.01) # compute the degeneracy
print("Degeneracy",deg)
print("Excited states",sc.get_excited(n=int(deg)+1,mode=mode))











