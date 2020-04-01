# Add the root path of the dmrgpy library and uncomment the next two lines
# for example PATH = /home/jose/programs/dmrgpy/src
#PATH = PATH_TO_DMRGPY_LIBRARY
#import os ; import sys ; sys.path.append(PATH)

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import spinchain
n = 4 # number of spin sites
import time
# first let us create a Hubabrd model

fc = fermionchain.Spinful_Fermionic_Chain(n) # create the object
#fc.maxm = 20
#fc.nsweeps = 40
U = 10.0

h0 = 0
for i in range(n-1):
    h0 = h0 + fc.Cdagup[i]*fc.Cup[i+1]
    h0 = h0 + fc.Cdagdn[i]*fc.Cdn[i+1]
h0 = h0 + h0.get_dagger()
for i in range(n):
   # h0 = h0 + U*fc.Nup[i]*fc.Ndn[i] -U*0.5*(fc.Nup[i]+fc.Ndn[i])
    h0 = h0 + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-0.5)

fc.set_hamiltonian(h0)
print([fc.vev(o).real for o in fc.Sx])
print([fc.vev(o).real for o in fc.Sy])
print([fc.vev(o).real for o in fc.Sz])
#print(fc.get_excited(n=4,mode="ED"))
import os

es = np.linspace(-1,6,600)

t0 = time.time()
ii = 1
jj = 1
name = (fc.Sx[ii],fc.Sx[jj])
(esh,ch) = fc.get_dynamical_correlator(name=name,mode="DMRG",es=es)
scale = 4./U
esh /= scale # scale by the estimated effective exchange
ch *= scale
t1 = time.time()
print("Time with fermions",t1-t0)
# now create the Heisenberg chain
sc = spinchain.Spin_Chain([2 for i in range(n)])
h1 = 0
for i in range(n-1):
    h1 = h1 + sc.Sx[i]*sc.Sx[i+1]
    h1 = h1 + sc.Sy[i]*sc.Sy[i+1]
    h1 = h1 + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h1)
name = (sc.Sx[ii],sc.Sx[jj])
(ess,cs) = sc.get_dynamical_correlator(name=name,mode="DMRG",es=es)
t2 = time.time()
print("Time with spins",t2-t1)

import matplotlib.pyplot as plt

plt.scatter(esh,ch.real,c="red",label="Hubbard")
plt.plot(ess,cs.real,c="blue",label="Heisenberg")
#plt.plot(inds,mz2,c="green",label="Exact")
plt.legend()
plt.xlabel("Site")
plt.ylabel("Correlator")
plt.show()



