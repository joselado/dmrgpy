# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create an S=1/2 spin chain
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange (Heisenberg model)
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

# you may need to make this bigger for bigger chains
sc.kpmmaxm = 20 # band dimension in KPM
sc.maxm = 20 # bond dimension in GS

xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
xs = range(n) # indexes
ys = np.linspace(-0.5,3.0,200) # frequencies
def compute(i):
  print("Computing ",i)
  name = (sc.Sz[i],sc.Sz[0]) # these are the two operators of the dyn. corr.
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",name=name,
          es=ys,delta=0.05)
  return s

zs = [compute(x) for x in xs] # compute all the dynamical correlators

zs = np.array(zs).T # transpose the array

# plot the results
import matplotlib.pyplot as plt
import matplotlib

# the Fourier transform in the x direction would be S(q,w)
matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.contourf(xs,ys,zs,100,cmap="rainbow")
plt.ylabel("Frequency [J]")
plt.xlabel("Distance")
plt.colorbar(label="S(r,w)")
plt.show()












