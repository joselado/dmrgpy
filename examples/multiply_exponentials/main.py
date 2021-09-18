# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain
#sc.itensor_version = "julia"

hfe = sum(sc.Sz) # Fully ferromagnetic
sc.set_hamiltonian(-hfe) ; wf = sc.get_gs() # Fully up wavefunction

wf0 = wf.copy() # copy the initial wavefunction

# rotate around the x axis 180 degrees (so that every site becomes orthogonal)
for i in range(n): 
    wf = sc.exponential(1j*sc.Sx[i]*np.pi,wf)

print("Overlap after rotation",wf0.overlap(wf))










