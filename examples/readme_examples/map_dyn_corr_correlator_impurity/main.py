# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(4)] # spin 1/2 heisenberg chain
spins = spins + ["S=1"] + spins # put S=1 in the middle
n = len(spins) # total number of spins
# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
for i in range(n): # loop over sites
  name = (sc.Sz[i],sc.Sz[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s.real) # store

xs = range(n)
ys = e
zs = np.array(zs).transpose()

# plot the results
import matplotlib.pyplot as plt
import matplotlib
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})


matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
zs = zs/np.max(zs)
plt.contourf(xs,ys,zs,100,cmap="magma",vmax=0.8)
plt.ylabel("frequency [J]")
plt.xlabel("Site")
plt.show()












