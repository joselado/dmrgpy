# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 3
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)]) # create the chain
h = 0


for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
h = h + sc.Sx[0]*sc.Sx[n-1]
h = h + sc.Sy[0]*sc.Sy[n-1]
h = h + sc.Sz[0]*sc.Sz[n-1]

h = h + h.get_dagger()

sc.set_hamiltonian(h)

sc.get_gs(mode="ED") # compute ground state

### Compute the dyamical correlator ###

i,j = np.random.randint(n),np.random.randint(n)
j = i
name = (sc.Sy[i],sc.Sy[j]) # operator

es = np.linspace(-0.5,6.0,100) # energies of the correlator
delta = 1e-1 # smearing of the correlator

#submodes = ["CVM","INV","ED","EX","KPM"] # all ED submodes
submodes = ["CVM","INV","ED"] # selected ED submodes

s = len(submodes)*4 # sizes
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(6,3))
import time
for submode in submodes: # loop over submodes
    t0 = time.time()
#    x,y = sc.get_dynamical_correlator(mode="ED",name=name,submode=submode,
#        es=es,delta=delta)
    from dmrgpy.dynamicstk import spincorrelators
    x,y = spincorrelators.get_full_SS_correlator(sc,mode="ED",
            submode=submode,i=i,es=es,delta=1e-2)
    t1 = time.time()
    label = submode + ", T="+str(np.round(t1-t0,1)) # label
    plt.plot(x,y.real,label=label,markersize=s,marker="o")
    s = s*0.5

plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.tight_layout()
plt.show()










