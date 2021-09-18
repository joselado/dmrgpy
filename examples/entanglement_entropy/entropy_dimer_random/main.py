# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain

from dmrgpy import mps

for i in range(10):
    wf = (mps.random_mps(sc) + mps.random_mps(sc)).normalize()
    print(wf.get_site_entropy(0))









