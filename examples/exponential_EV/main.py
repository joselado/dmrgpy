# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

Mz = sum(sc.Sz) # total magnetization in Z
Mx = sum(sc.Sx) # total magnetization in X

def get(mode):
    sc.set_hamiltonian(Mz) # only magnetic field
    wf = sc.get_gs(mode=mode) # Z ferromagnetic wavefunction
    # The operator in the exponent MUST be Hermitian
    wf1 = sc.exponential(Mx,wf,mode=mode) # compute exp(Mx)*wf
    c = sc.overlap(wf,wf1,mode=mode) # compute <wf|exp(Mx)|wf>
    return c # return the number

# Compare the result between DMRG and ED
print("ED",get("ED"))
print("DMRG",get("DMRG"))
