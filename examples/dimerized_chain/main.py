# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
import matplotlib.pyplot as plt
djs = np.linspace(0.0,0.1,40)
fo = open("SWEEP.OUT","w")
for dj in djs:
  print("Doing",dj)
  n = 60
  spins = [2 for i in range(n)]
  sc = spinchain.Spin_Hamiltonian(spins) # create the chain
  
  b = 1.0+dj
  phi = np.pi*np.random.random()
  
  def fj(i,j):
    if i==(j-1):
      return (1.0 + dj*(-1)**i)*np.identity(3)
#      return (1.0 + 0.3*np.cos(i*b*np.pi + phi))*np.identity(3)
    else: return 0.0
  
  sc.set_exchange(fj) # create couplings
  
  def fm(i):
#    return dj*np.cos(i*b)*np.array([0.,0.,1.])
    if i==n//3: return [3.0,0.,0.]
    else: return [0.,0.,0.]
  sc.set_fields(fm)
  pairs = [(n//3,n//3 + i) for i in range(1,2*n//3)]
#  cs = sc.correlator(pairs) # get the correlator
  cs = sc.magnetization()[0] # compute the magnetization
  cs = cs[n//3:2*n//3]
  
  plt.plot(range(len(cs)),cs,label=str(dj)) # correlator using DMRG
  for i in range(len(cs)):
    fo.write(str(dj)+"   ")
    fo.write(str(i)+"   ")
    fo.write(str(abs(cs[i]))+"\n")
  fo.flush()
fo.close()
plt.legend()
plt.show()


