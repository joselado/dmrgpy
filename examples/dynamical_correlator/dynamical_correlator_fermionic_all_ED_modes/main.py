# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0


for i in range(n):
  for j in range(n):
      t = np.random.random() + 1j*np.random.random()
      h = h + fc.Cdag[i]*fc.C[j]*t # random hopping
      V = np.random.random()*0.2
      h = h + V*fc.N[i]*fc.N[j] # random interaction


h = h + h.get_dagger()

fc.set_hamiltonian(h)

fc.get_gs() # compute ground state

### Compute the dyamical correlator ###

i,j = np.random.randint(n),np.random.randint(n)
op = fc.Cdag[i]*fc.C[j] # operator
name = (op.get_dagger(),op) # operator

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
    x,y = fc.get_dynamical_correlator(mode="ED",name=name,submode=submode,
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










