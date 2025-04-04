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
    ents = [] # entanglement entropies
    ts = [] # times
    for i in range(nt): # loop over time steps
        ents.append(wf.get_correlation_entropy()) # compute correlation entropy
        ts.append(i*dt) # store time
        wf = fc.exponential(1j*dt*h,wf) # apply e^{iHdt}|WF>
    return ts,ents


import matplotlib.pyplot as plt

Vs = [0.,1.,2.] # interactions

for V in Vs:
    (ts,ents) = get_entropy(V=V)
    plt.plot(ts,ents,marker="o",label="V = "+str(V))

plt.xlabel("Time")
plt.ylabel("Correlation entropy")
plt.ylim([-0.1,1.])
plt.legend()
plt.show()











