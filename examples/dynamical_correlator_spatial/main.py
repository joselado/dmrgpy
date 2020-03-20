# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 8
# create an S=1/2 spin chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
def fj(i,j): # first neighbor coupling
  if abs(i-j)==1: 
#      return np.array([[1,0,0],[0,1,0],[0,0,0]])
      return 1.0
  else: return 0.0


sc.set_exchange(fj)

sc.kpmmaxm = 20 # KPM maxm
xs = [] # empty list
ys = [] # empty list
zs = [] # empty list
for i in range(n): # loop over sites
  print("Doing ",i)
  (e,s) = sc.get_dynamical_correlator(mode="DMRG",i=i,j=i,name="XX",
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



