# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

hfe = sum(sc.Sz) # Fully ferromagnetic
sc.set_hamiltonian(-hfe) ; wfup = sc.get_gs() # Fully up wavefunction
sc.set_hamiltonian(hfe) ; wfdn = sc.get_gs() # Fully down wavefunction

print("Overlap between orthogonal wavefunctions",wfup.dot(wfdn).real)

# Flip all the spin in each site, to go from up to down
wfflip = wfdn.copy() 
for i in range(n): wfflip = 2.*sc.Sx[i]*wfflip

# Now do some simple algebra with MPS
print("Rotate and overlap",wfup.dot(wfflip).real)
print("Magnetization up",wfup.dot(hfe*wfup).real)
print("Magnetization down",wfdn.dot(hfe*wfdn).real)

wfs = wfup + 3.7*wfdn # sum the two wavefunctions

print("Overlap with the sum",wfdn.overlap(wfs).real)










