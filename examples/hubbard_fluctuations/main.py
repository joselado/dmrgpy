# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 10 # number of spinful fermionic sites
def get_fc(u):
  fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
  h = 0
  for i in range(n-1):
      h = h + fc.Cdagup[i]*fc.Cup[i+1]
      h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
  for i in range(n):
      h = h + u*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
  h = h + h.get_dagger()
  fc.set_hamiltonian(h) 
  
  # Compute the dynamical correlator defined by
  # <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>
  
  fc.nsweeps = 6
  fc.maxm = 20 # maximum bond dimension in KPM
  return fc

us = np.linspace(-10.0,10.0,20)
ds = []
for u in us:
  fc = get_fc(u)
  d = np.mean(fc.get_density_fluctuation())
  ds.append(d)
  print(u,d)
ds = np.array(ds)
np.savetxt("D_VS_MU.OUT",np.matrix([us,ds]).T)
import matplotlib.pyplot as plt
plt.plot(us,ds,marker="o")
plt.show()
# The result will be written in a file called DYNAMICAL_CORRELATOR.OUT
fc.gs_energy()
exit()











