# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the chain
m = np.random.random((n,n)) + 1j*np.random.random((n,n))

m = m + np.conjugate(m).T

h = 0

for i in range(n):
  for j in range(n):
      h = h + m[i,j]*fc.Cdag[i]*fc.C[j]


fc.set_hamiltonian(h) # Hamiltonian
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization


### Compute the dyamical correlator ###

ii,jj = 1,1
ni,nj = fc.vev(fc.N[ii]),fc.vev(fc.N[jj])
name = (fc.N[ii]-ni*fc.Id,fc.N[jj]-nj*fc.Id)

es = np.linspace(-0.5,6.0,500) # energies of the correlator
delta = 1e-1 # smearing of the correlator
import time
t1 = time.time()
x1,y1 = fc.get_dynamical_correlator(mode="DMRG",name=name,es=es,delta=delta)

import scipy.linalg as lg

(ei,vi) = lg.eigh(m)
vi = vi.T # transpose
eta = 6e-2
y2 = y1*0.0
for i in range(n):
  for j in range(n):
      if ei[i]<0.0 and ei[j]>0.0:
        v = np.abs(vi[i][ii])*np.abs(vi[j][jj])
        dw = es+ei[i]-ei[j]
        y2 +=  eta/(eta**2 + dw**2)*v**2/np.pi


### Plot the result ###

import matplotlib.pyplot as plt

#plt.plot(x0,y0.real,label="ED",marker="o")
plt.plot(x1,y1.real,label="DMRG",marker="o")
plt.plot(x1,y2.real,label="Lindhard",marker="o")
#plt.plot(x2,y2.real,label="DMRG-TD",marker="o")
#plt.plot(x3,y3.real,label="ED",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()










