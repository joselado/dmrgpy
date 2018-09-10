import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg

import spinchain

ns = np.array(range(4,30,2))
es = []

n = 60
spins = [3 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain

m = sc.sisj_edge()
  
import matplotlib.pyplot as plt

plt.plot(m[0],m[1],marker="o")

plt.show()
