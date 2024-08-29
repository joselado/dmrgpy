# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4
sc = fermionchain.Spinful_Fermionic_Chain(n) # create the spin chain
delta = -0.2 # dimerization
ks = np.linspace(0.,1.,20,endpoint=False)
wfs = []

for k in ks: # loop over kpoints
    h = 0
    for i in range(n-1):
        Ji = 1. + delta*(-1)**i # exchange
        h = h + Ji*sc.Cdagup[i]*sc.Cup[i+1]
        h = h + Ji*sc.Cdagdn[i]*sc.Cdn[i+1]
    h = h.get_dagger()
    z = np.exp(1j*np.pi*k) # complex twist
    zc = np.conjugate(z) # conjugate
    # add the last link
    Ji = 1. + delta*(-1)**(n-1) # exchange
    hk = z*Ji*(sc.Cup[0]*sc.Cup[n-1])
    hk = hk + zc*Ji*(sc.Cdn[0]*sc.Cdn[n-1])
    h = h + hk + hk.get_dagger()
    # set the Hamiltonian
    sc.set_hamiltonian(h)
    wf0 = sc.get_gs(mode="ED") # compute ground state
    wfs.append(wf0) # store

phi = 1. # initialize
for i in range(len(wfs)-1): # all pairs
    phi = phi*wfs[i].dot(wfs[i+1])
phi = phi*wfs[len(wfs)-1].dot(wfs[0]) # last link

print("Zak phase",np.angle(phi)/np.pi)







