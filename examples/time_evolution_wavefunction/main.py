# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 2 # number of sites in the chain
spins = ["S=1/2" for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins) # create the chain

# first Hamiltonian, to create an initial state
h0 = 0
for i in range(n): h0 = h0+(-1)**i*sc.Sz[i]
sc.set_hamiltonian(h0) # set AF hamiltonian
wf0 = sc.get_gs() # compute ground state

# second Hamiltonian, this will make the time evolution
h = 0
for i in range(n-1):
    h = h+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1] 
h = (h +h.get_dagger())


from dmrgpy.timeevolution import imaginary_exponential as evolve

def evolve_and_entropy(t):
    """Evolve for some time and compute the entropy of the resulting state"""
    wf = evolve(h,wf0,ts=[t])[0] # this is the final wavefunction at t
    return wf.get_site_entropy(0) # comptue the entropy

# compute entropies for different times
ts = np.linspace(0.,4.,40)
ss = []
for t in ts:
    s = evolve_and_entropy(t)
    print(t,s)
    ss.append(s)


import matplotlib.pyplot as plt


plt.plot(ts,ss,marker="o")
plt.xlabel("Time")
plt.ylabel("Entropy")
plt.show()
