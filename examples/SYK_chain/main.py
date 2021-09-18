# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 4
L = 3
fc = fermionchain.Fermionic_Chain(n*L) # create the chain
h = 0
#for i in range(n):
#  for j in range(i):
#    h = h + np.random.random()*fc.C[i]*fc.Cdag[j]
C = fc.C
Cdag = fc.Cdag
for ll in range(L):
  for i in range(n):
    for j in range(i+1,n):
      for k in range(n):
        for l in range(k+1,n):
          Ni = C[i+ll*n]*C[j+ll*n]
          Nj = Cdag[k+ll*n]*Cdag[l+ll*n]
          h = h + (np.random.random()-.5)*Ni*Nj
          if ll<(L-1):
            Ni = C[i+ll*n]*C[j+ll*n]
            Nj = Cdag[k+(ll+1)*n]*Cdag[l+(ll+1)*n]
            h = h + (np.random.random()-.5)*Ni*Nj
h = h + h.get_dagger()

ms = [3,10,20,1000]
#ms = [10]
#fc.set_hamiltonian(h)
#e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
#print(e0)
import profile
#
#def f():
#  fc.set_hamiltonian(h)
#  e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
#
#profile.run('f()','restats')
#import pstats
#from pstats import SortKey
#p = pstats.Stats('restats')
#p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
#exit()
import time
fc.set_hamiltonian(h)
t0 = time.time()
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
t1 = time.time()
print("Time in GS = ",t1-t0)
print(np.mean(fc.get_density()))
exit()
#fc.mpomaxm = 500
wf = fc.get_gs()
#wf = fc.C[0]*wf
#fc.applyinverse(h-e1-1j*1e-1,wf)
(x,y) = fc.get_dynamical_correlator(name=(fc.Cdag[0],fc.C[0]),submode="TD")
#t2 = time.time()
print("Time in KPM = ",t2-t1)
exit()
import matplotlib.pyplot as plt
plt.plot(x,y)
plt.show()
exit()
t2 = time.time()
print("Time in correlator = ",t2-t1)










