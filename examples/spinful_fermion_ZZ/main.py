# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import spinchain
n = 6 # number of spin sites

# first let us create a Hubabrd model

fc = fermionchain.Spinful_Fermionic_Chain(n) # create the object
fc.maxm = 20
U = 1.0
def ft(i,j):
    if abs(i//2-j//2)==1 and i%2==j%2: return 1.0 # first neighbor coupling
    if i==j: return -2*U # set to half filling
    return 0.0

def fu(i,j):
    if i==j: return U
    return 0.0

def fxc(i,j):
    if abs(i-j)==1: return np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,1.]])
    else: return 0.0


fc.set_hoppings(ft) # add the term to the Hamiltonian
fc.set_hubbard(fu) # add the term to the Hamiltonian
fc.set_exchange(fxc) # add the term to the Hamiltonian
pairs = [(0,i) for i in range(n)]
for mode in ["ED","DMRG"]:
  print("Mode ",mode)
  ch = fc.get_correlator(pairs=pairs,name="ZZ",mode=mode).real
  print("Total density")
  print(fc.get_density_spinful(mode=mode))
  print("YY correlator")
  print(ch)
  print("ground state energy")
  print(fc.gs_energy(mode=mode))
