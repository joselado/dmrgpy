import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg


import spinchain

ns = np.array(range(4,30,2))
es = []

for n in ns:
  spins = [2 for i in range(n)]
  sc = spinchain.Spin_Hamiltonian(spins) # create the chain
  e = sc.gs_energy()
  es.append(e)
  print(e)

es = np.array(es) # array
  
import matplotlib.pyplot as plt

plt.scatter(ns,es/ns)

plt.show()
