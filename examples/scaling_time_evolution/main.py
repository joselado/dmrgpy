# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
import time
ns = np.array(range(4,30,2))
es = []
n = 4
def compute(n):
  spins = [2 for i in range(n)]
  sc = spinchain.Spin_Hamiltonian(spins) # create the chain
  sc.clean()
  sc = spinchain.Spin_Hamiltonian(spins) # create the chain
  def fj(i,j):
      if 0.9<abs(i-j)<1.1: return 1.0
      return 0.0
  sc.set_exchange(fj) # set exchange couplings
  #sc.set_fields(lambda x: [0.2,0.2,0.2]) # set exchange couplings
  sc.maxm = 10
  sc.nsweeps = 2
  sc.get_gs()
  i = 0 ; j=1
  t0 = time.time()
  sc.fit_td = True
#  (x,y) = sc.evolution(nt=100,dt=0.03,name="ZZ",i=i,j=j)
  (x2,y2) = sc.get_dynamical_correlator(es=es,use_kpm=False,dt=0.1)
  t1 = time.time()
  print("Time",t1-t0)
  return t1-t0
ns = np.array(range(10,12,2))
ns = [6]
ts = np.array([compute(n) for n in ns])
out = np.polyfit(np.log(ns),np.log(ts),deg=1)
print("Power",out[0])
import matplotlib.pyplot as plt
plt.plot(np.log(ns),np.log(ts),marker="o")
plt.ylabel("Log Time")
plt.xlabel("Log Size")
#plt.plot(x,y.imag)
#plt.scatter(x,y.imag)
plt.show()
