# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
ns = np.array(range(4,30,2))
es = []
n = 8
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain

h = 0 
for i in range(n-1): h = h +fc.Cdag[i]*fc.C[i+1]
#for i in range(n-1): h = h +fc.Cdag[i]*fc.Cdag[i+1]
for i in range(n-1): h = h +4*fc.N[i]*fc.N[i+1]
h = h + h.get_dagger()
fc.set_hamiltonian(h)


A = fc.Cdag[0]
B = fc.N[0]
nt = 1e3 # number of time steps
dt = 1e-2 # time step

from dmrgpy import timedependent

(x,y) = timedependent.evolution_ABA(fc,nt=nt,dt=dt,mode="ED",A=A,B=B)
(x1,y1) = timedependent.evolution_ABA(fc,nt=nt,dt=dt,mode="DMRG",A=A,B=B)
import matplotlib.pyplot as plt
plt.subplot(1,2,1)
plt.plot(x,y.real,c="red",label="ED")
plt.scatter(x1,y1.real,c="blue",label="DMRG")
plt.subplot(1,2,2)
plt.plot(x,y.imag,c="red",label="ED")
plt.scatter(x1,y1.imag,c="blue",label="DMRG")
plt.legend()
#plt.plot(x,y.imag)
#plt.scatter(x,y.imag)
plt.show()


