# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain


def get_entropy(V=0):
    n = 4
    fc = fermionchain.Fermionic_Chain(n) # create the chain
    h0,hv = 0 ,0
    for i in range(n-1): 
        h0 = h0 +fc.Cdag[i]*fc.C[i+1] # hopping
    for i in range(n-1): 
        hv = hv +(fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    h0 = h0 + h0.get_dagger()
    h = h0 + V*hv
    fc.set_hamiltonian(h)
    
    nt = 100 # number of time steps
    dt = 1e-1 # time step
    
    mode="ED"
    wf0 = fc.get_gs(mode=mode) # compute ground state
    wf = wf0.copy() # copy wavefunction
    wf = fc.Cdag[1]*fc.C[0]*wf0 # create an electron-hole pair
    wf = wf.normalize() # normalize wavefunction
    ents = [] # entanglement entropies
    occs0 = []
    occs1 = []
    ts = [] # times
    for i in range(nt): # loop over time steps
        ents.append(wf.get_correlation_entropy()) # compute correlation entropy
        occs0.append(wf.dot(fc.N[0]*wf)) # occupation in site 0
        occs1.append(wf.dot(fc.N[1]*wf)) # occupation in site 1
        ts.append(i*dt) # store time
        wf = fc.exponential(1j*dt*h,wf) # apply e^{iHdt}|WF>
    return ts,ents,occs0,occs1


import matplotlib.pyplot as plt

Vs = [0.,1.,2.] # interactions

fig = plt.figure(figsize=(12,4))

for V in Vs:
    (ts,ents,occs0,occs1) = get_entropy(V=V)
    plt.subplot(1,3,1)
    plt.plot(ts,ents,marker="o",label="V = "+str(V))
    plt.subplot(1,3,2)
    plt.plot(ts,occs0,marker="o",label="V = "+str(V))
    plt.subplot(1,3,3)
    plt.plot(ts,occs1,marker="o",label="V = "+str(V))
plt.subplot(1,3,1)
plt.xlabel("Time"); plt.ylabel("Correlation entropy") ; plt.ylim([-0.1,1.])
plt.legend()

plt.subplot(1,3,2)
plt.xlabel("Time"); plt.ylabel("Occupation in 0") ; plt.ylim([-0.1,1.1])
plt.legend()

plt.subplot(1,3,3)
plt.xlabel("Time"); plt.ylabel("Occupation in 1") ; plt.ylim([-0.1,1.1])
plt.legend()


plt.tight_layout()
plt.show()











