# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
import meanfield
ns = np.array(range(4,30,2))
es = []
n = 6
spins = [2 for i in range(n)]
fo = open("SWEEP.OUT","w")
ps = np.linspace(0.,1.,20)
for p in ps: # loop over poarameters
  sc = spinchain.Spin_Chain(spins) # create the chain
  m0 = [[0.,0.,(-1)**i] for i in range(sc.ns)]
  sc = meanfield.spinchain_meanfield(sc,m0=m0,p=p) # return spin chain
  (mx,my,mz) = sc.get_magnetization() # get magnetization
  for i in range(1): # loop
      fo.write(str(p)+"   ")
#      fo.write(str(i)+"   ")
      fo.write(str(mz[i])+"\n")
  fo.flush()
fo.close()











