# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 8
# create an S=1/2 spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)


sc.kpmmaxm = 20 # KPM maxm
xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
for i in range(n): # loop over sites
  print("Doing ",i)
  name = (sc.Sx[i],sc.Sx[i])
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=np.linspace(-0.5,4.0,200),delta=0.05)
  zs.append(s.real) # store

xs = range(n)
ys = e
zs = np.array(zs).transpose()

# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.contourf(xs,ys,zs,100)
plt.ylabel("frequency [J]")
plt.xlabel("Site")
plt.show()



