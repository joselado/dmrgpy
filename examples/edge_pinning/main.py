# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,30,2))
es = []
n = 20
spins = [2 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
def fm(i):
  if i==0: return [3.0,0.,0.]
  else: return [0.,0.,0.]
sc.set_fields(fm) # add those local exchange fields
cs = sc.get_magnetization()[0] # compute the magnetization
  
import matplotlib.pyplot as plt
plt.plot(range(len(cs)),cs) # correlator using DMRG
plt.show()


