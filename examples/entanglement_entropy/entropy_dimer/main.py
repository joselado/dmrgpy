# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain

# Heisenberg model
h = 0
h = h +sc.Sx[0]*sc.Sx[1]
h = h +sc.Sy[0]*sc.Sy[1]
h = h +sc.Sz[0]*sc.Sz[1]

# Add a magnetic field in the first site to polarize it
bzs = np.linspace(0.,4.0,10) # magnetic field
ss = [] # lists for entropies
xx,yy,zz = [],[],[] # lists for correlators

Oxx = sc.Sx[0]*sc.Sx[1] # define operator
Oyy = sc.Sy[0]*sc.Sy[1] # define operator
Ozz = sc.Sz[0]*sc.Sz[1] # define operator


for bz in bzs:
    h1 = h + bz*sc.Sz[0] # add magnetic field
    sc.set_hamiltonian(h1)
    wf = sc.get_gs() # compute ground state
    s = wf.get_site_entropy(0) # compute entanglement entropy of the first site
    ss.append(s) # store
    xx.append(wf.dot(Oxx*wf).real) # store correlator
    yy.append(wf.dot(Oyy*wf).real) # store correlator
    zz.append(wf.dot(Ozz*wf).real) # store correlator

import matplotlib.pyplot as plt

plt.subplot(121) # plot entropy
plt.plot(bzs,ss,marker="o")
plt.xlabel("Local Bz in site #1")
plt.ylabel("Entanglement entropy")


plt.subplot(122) # plot correlators
plt.plot(bzs,xx,marker="o",label="$\\langle S^x_1 S^x_2 \\rangle$")
plt.plot(bzs,yy,marker="o",label="$\\langle S^y_1 S^y_2 \\rangle$")
plt.plot(bzs,zz,marker="o",label="$\\langle S^z_1 S^z_2 \\rangle$")
plt.xlabel("Local Bz in site #1")
plt.ylabel("Spin correlator")
plt.legend()
plt.show()










