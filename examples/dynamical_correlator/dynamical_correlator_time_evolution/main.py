# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,30,2))
es = []
n = 5
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

sc.set_hamiltonian(h)

sc.maxm = 10
es = np.linspace(-1.0,10.0,4000)
name = (sc.Sx[0],sc.Sx[0]) # correlator to compute


import time
print("Starting")
t0 = time.time()
(x1,y1) = sc.get_dynamical_correlator(es=es,name=name,submode="KPM")
t1 = time.time()
sc.fit_td = True
sc.tevol_custom_exp = True
(x2,y2) = sc.get_dynamical_correlator(es=es,name=name,submode="TD")
t2 = time.time()
print("Time in TD",t2-t1)
print("Time in KPM",t1-t0)
import matplotlib.pyplot as plt
plt.plot(x1,y1.real,label="KPM")
plt.plot(x2,np.abs(y2),label="TD")
plt.legend()
#plt.scatter(x,y.imag)
plt.show()











