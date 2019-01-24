# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,30,2))
es = []
n = 20
spins = [3 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
m = sc.sisj_edge()
  
import matplotlib.pyplot as plt
plt.plot(m[0],m[1],marker="o")
plt.show()


