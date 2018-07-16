from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 20
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def fb(i): return [0,0,np.cos(i*np.pi*np.sqrt(2))*(-1)**i]
sc.set_fields(fb)

#e = sc.get_excited(n=10)
#print(e)
#exit()



fo = open("DCF.OUT","w") # dynamical correlation function

for i in range(n): # loop over sites
  (xs,ys) = sc.get_spismj(n=2000,mode="DMRG",i=i,j=i)
  print("Doing",i)
  for (x,y) in zip(xs,ys):
    fo.write(str(i)+"  ")
    fo.write(str(x)+"  ")
    fo.write(str(y)+"\n")
  fo.flush()

fo.close()



