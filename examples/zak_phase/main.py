# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 4
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
delta = -0.2 # dimerization
ks = np.linspace(0.,1.,20,endpoint=False)
wfs = []

for k in ks: # loop over kpoints
    h = 0
    for i in range(n-1):
        Ji = 1. + delta*(-1)**i # exchange
        h = h + Ji*sc.Sx[i]*sc.Sx[i+1]
        h = h + Ji*sc.Sy[i]*sc.Sy[i+1]
        h = h + Ji*sc.Sz[i]*sc.Sz[i+1]
    z = np.exp(1j*np.pi*2.*k) # complex twist
    zc = np.conjugate(z) # conjugate
    Spn = sc.Sx[n-1] + 1j*sc.Sy[n-1] # S+_n
    Sp0 = sc.Sx[0] + 1j*sc.Sy[0] # S+_n
    Smn = Spn.get_dagger() # S-_n
    Sm0 = Sp0.get_dagger() # S-_n
    # add the last link
    Ji = 1. + delta*(-1)**(n-1) # exchange
    h = h + Ji*(z*(Spn*Sm0)/2. + zc*(Sp0*Smn)/2. + sc.Sz[0]*sc.Sz[n-1])
    # set the Hamiltonian
    sc.set_hamiltonian(h)
    wf0 = sc.get_gs(mode="DMRG") # compute ground state
    wfs.append(wf0) # store

phi = 1. # initialize
for i in range(len(wfs)-1): # all pairs
    phi = phi*wfs[i].dot(wfs[i+1])
phi = phi*wfs[len(wfs)-1].dot(wfs[0]) # last link

print("Zak phase",np.angle(phi)/np.pi)







