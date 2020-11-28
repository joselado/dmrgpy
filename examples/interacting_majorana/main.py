# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 40 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
V = -0.9
for i in range(n-1): # hopping
    h = h + fc.Cdag[i]*fc.C[i+1]
#    h = h + 0.5*(-1)**i*fc.Cdag[i]*fc.C[i+1]
#    h = h + 0.3*fc.C[i]*fc.C[i+1]
    h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 40
fc.nsweeps = 10
fc.set_hamiltonian(h) # set the hoppings
fc.get_gs()

#print(fc.get_correlation_entropy()) ; exit()
out = fc.get_correlated_orbitals(ordered=False)
f = open("MAP.OUT","w")
for i in range(len(out)):
  for j in range(len(out[0])):
      f.write(str(i)+"  ")
      f.write(str(j)+"  ")
      f.write(str(out[i][j])+"\n")

f.close()
