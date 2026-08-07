# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain

from dmrgpy import mps

ss = []
for i in range(10):
    wf = (mps.random_mps(sc) + mps.random_mps(sc)).normalize()
    s = wf.get_site_entropy(0)
    print(s)
    ss.append(s)

import matplotlib.pyplot as plt

plt.plot(range(len(ss)),ss,marker="o")
plt.xlabel("random trial")
plt.ylabel("Site(0) entanglement entropy")
plt.title("Entanglement entropy of random 2-site MPS superpositions")
plt.show()









