# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n  = 2 # total number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain

tij = [[.5,1.0],[1.0,-0.5]] # hopping
H = 0
for i in range(n):
    for j in range(n):
        H = H +tij[i][j]*fc.Cdagup[i]*fc.Cup[j]
        H = H +tij[i][j]*fc.Cdagdn[i]*fc.Cdn[j]

fc.set_hamiltonian(H) # set the Hamiltonian


mode = "ED"
submode = "INV"
es = np.linspace(0.,4.,500)
wf = fc.get_gs(mode=mode)
wf1 = fc.Cdagup[0]*wf
wf2 = fc.Cdagup[1]*wf
print(wf1.dot(wf1))
print(wf2.dot(wf2))
T01 = (fc.Cup[0],fc.Cdagup[1])
T10 = (fc.Cup[1],fc.Cdagup[0])
A01 = fc.Cdagup[0]*fc.Cup[1]
A10 = fc.Cdagup[1]*fc.Cup[0]
A10 = A10 - wf.dot(A10*wf) # reddefine
A01 = A01 - wf.dot(A01*wf) # reddefine
T01 = (A01,A01.get_dagger())
T10 = (A10,A10.get_dagger())

es01,ds01 = fc.get_dynamical_correlator(name=T01,mode=mode,submode=submode,es=es)
es10,ds10 = fc.get_dynamical_correlator(name=T10,mode=mode,submode=submode,es=es)

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,4))
plt.plot(es01,ds01.real,label="0 -> 1",c="blue")
plt.plot(es10,ds10.real,label = "1 -> 0",c="red")
plt.xlabel("$\\omega$")
plt.ylabel("$A(\\omega)$")

plt.tight_layout()
plt.legend()
plt.show()


