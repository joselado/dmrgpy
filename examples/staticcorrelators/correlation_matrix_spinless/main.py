# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 60 # number of fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
U = -0.95
for i in range(n-1): # hopping
    h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 20
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

# out[i][j] is the weight of orbital i on site j: a natural heatmap of
# the correlated-orbital decomposition of the ground state
import matplotlib.pyplot as plt
plt.imshow(np.array(out),aspect="auto",origin="lower")
plt.colorbar(label="Orbital weight")
plt.xlabel("Site j")
plt.ylabel("Orbital i")
plt.show()





