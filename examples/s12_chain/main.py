# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,30,2))
es = []
for n in ns:
  spins = [2 for i in range(n)]
  sc = spinchain.Spin_Chain(spins) # create the chain
  e = sc.gs_energy()
  es.append(e)
  print(e)
es = np.array(es) # array
  
import matplotlib.pyplot as plt
plt.scatter(ns,es/ns)
plt.show()











