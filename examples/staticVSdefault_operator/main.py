# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 60 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

sc.set_hamiltonian(sum(sc.Sz)) ; wf = sc.get_gs() # Fully up wavefunction

h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]


Op = h # operator
from dmrgpy.multioperatortk.staticoperator import StaticOperator
SOp = StaticOperator(Op,sc) # create a static operator

import time
for (O,mode) in [(SOp,"static"),(Op,"default")]:
    t0 = time.time()
    for i in range(50):
        wf = O*wf
    t1 = time.time()
    print("Time in",mode,t1-t0)
