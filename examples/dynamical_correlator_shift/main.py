# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
import time

def get(a):
  n = 10
  # create a random spin chain
  spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
  
  
  # create first neighbor exchange
  sc = spinchain.Spin_Chain(spins) # create the spin chain
  #sc.itensor_version = "julia"
  h = 0
  for i in range(n-1):
      h = h + sc.Sx[i]*sc.Sx[i+1]
      h = h + sc.Sy[i]*sc.Sy[i+1]
      h = h + sc.Sz[i]*sc.Sz[i+1]
  sc.set_hamiltonian(h)
  
  sc.kpmmaxm = 20 # KPM maxm
  sc.maxm = 20 # KPM maxm
  A,B = sc.Sz[1],sc.Sz[1]
  es = np.linspace(-0.5,6,2000)
  delta = 1e-2
  #sc.kpm_extrapolate = False
  e = sc.gs_energy()
  sc.kpm_extrapolate = True
  sc.kpm_extrapolate_factor = 2.0
  t0 = time.time()
  (x,y) = sc.get_distribution(X=h-e,xs=es,delta=delta,A=A,B=B,a=a,b=0.9)
  t1 = time.time()
  print(a,t1-t0)
  return (x,y)




# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
for a in np.linspace(-0.9,0.5,5):
  (x,y) = get(a)
  plt.plot(x,y.real,label=str(a))
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



