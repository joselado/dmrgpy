# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain
n = 20 # number of fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0
U = -0.95
for i in range(n-1): # hopping
    h = h + 1j*fc.Cdag[i]*fc.C[i+1] # complex, just to check things are ok
    h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
h = h + h.get_dagger()
##############################
# Setup the Many Body Hamiltonian
fc.maxm = 20
fc.nsweeps = 10
fc.set_hamiltonian(h) # set the hoppings
wf = fc.get_gs(mode="DMRG")

#from dmrgpy.mps import random_mps
#wf = random_mps(fc)

import time

t0 = time.time()
m1 = wf.get_correlation_matrix(dmmode="fast")
t1 = time.time()


m2 = wf.get_correlation_matrix(dmmode="full")
t2 = time.time()

print("Python mode",t1-t0)
print("C++ mode",t2-t1)

print("Difference",np.mean(np.abs(m1-m2)))

# compare the "fast" (Python) and "full" (C++) correlation matrices, and
# their difference, as heatmaps
import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,3,figsize=(12,4))
for ax,m,title in zip(axes,[m1,m2,np.abs(m1-m2)],
        ["Python (fast)","C++ (full)","|Python - C++|"]):
    im = ax.imshow(np.abs(m),aspect="auto",origin="lower")
    ax.set_title(title)
    ax.set_xlabel("Site j")
    ax.set_ylabel("Site i")
    fig.colorbar(im,ax=ax)
plt.tight_layout()
plt.show()

