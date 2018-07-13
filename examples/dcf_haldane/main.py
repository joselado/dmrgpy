import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 30
spins = [3 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

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



