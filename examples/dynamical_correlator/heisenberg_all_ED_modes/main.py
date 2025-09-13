# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 4
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)]) # create the chain
h = 0


for i in range(n):
  for j in range(n):
      h = h + sc.Sx[i]*sc.Sx[j]*np.random.random() # random exchange
      h = h + sc.Sy[i]*sc.Sy[j]*np.random.random() # random exchange
      h = h + sc.Sz[i]*sc.Sz[j]*np.random.random() # random exchange


h = h + h.get_dagger()

sc.set_hamiltonian(h)

sc.get_gs(mode="ED") # compute ground state

### Compute the dyamical correlator ###

i,j = np.random.randint(n),np.random.randint(n)
name = (sc.Sy[i],sc.Sy[j]) # operator

es = np.linspace(-0.5,6.0,1000) # energies of the correlator
delta = 2e-2 # smearing of the correlator

#submodes = ["CVM","INV","ED","EX","KPM"] # all ED submodes
submodes = ["CVM","INV","ED"] # selected ED submodes

s = len(submodes)*4 # sizes
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,3))
import time
for submode in submodes: # loop over submodes
    t0 = time.time()
    x,y = sc.get_dynamical_correlator(mode="ED",name=name,submode=submode,
        es=es,delta=delta)
    t1 = time.time()
    label = submode + ", T="+str(np.round(t1-t0,1)) # label
    plt.plot(x,y.real,label=label,markersize=s,marker="o")
    s = s*0.5

plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.tight_layout()
plt.show()










