# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
n = 20
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
# dimerized coupling
def fj(i,j):
  if abs(i-j)==1: 
    ij = (i+j)%4
    dj = -0.2
    if ij==1: return 1.0 + dj
    else: return 1.0 - dj
  return 0.0
sc.set_exchange(fj) # set those exchange couplings
sc.kpmmaxm = 10 # KPM max m
fo = open("DCF.OUT","w") # dynamical correlation function
for i in range(n): # loop over sites
  (xs,ys) = sc.get_spismj(n=1000,mode="DMRG",i=i,j=i)
  print("Doing",i)
  for (x,y) in zip(xs,ys):
    fo.write(str(i)+"  ")
    fo.write(str(x)+"  ")
    fo.write(str(y)+"\n")
  fo.flush()
fo.close()
