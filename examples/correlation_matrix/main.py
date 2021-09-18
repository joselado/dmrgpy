# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 20 # number of spinful fermionic sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
h = 0
U = 2.0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1]
    h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
for i in [0]: # Hubbard
    h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
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
exit()

cm = fc.get_correlation_matrix()
import scipy.linalg as lg
print(lg.eigvalsh(cm))









