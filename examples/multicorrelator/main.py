# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 8
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def fj(i,j):
  if abs(i-j)==1:
#      if ((i+j)/2.-0.5)%2==0.0: return 0.4
      return 1.0
  else: return 0.0
sc.set_exchange(fj)
from dmrgpy import multicorrelator
es = np.linspace(0.0,4.0,100)
multicorrelator.multicorrelator(sc,es=es)


