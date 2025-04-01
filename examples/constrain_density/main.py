# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 5
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
m = np.random.random((n,n)) + 1j*np.random.random((n,n))

m = m + np.conjugate(m).T

h = 0

for i in range(n):
  for j in range(n):
      h = h + m[i,j]*fc.Cdag[i]*fc.C[j]


Ntarget = 5

lamb = 100 # lagrange multiplier for the density
Ntot = sum(fc.Nup) + sum(fc.Ndn) - Ntarget # operator
Ntot = sum(fc.N) - Ntarget # operator
Ntot2 = Ntot*Ntot


#m = fc.get_ED_obj().MO2matrix(Ntot2)
#import scipy.linalg as lg
#print(lg.eigvalsh(m.todense()))
#exit()


h = h - lamb*(Ntot*Ntot) # add the lagrange multiplier
h = lamb*Ntot2 # add the lagrange multiplier

fc.set_hamiltonian(h) # Hamiltonian
wf = fc.get_gs(mode="ED") # energy with exact diagonalization

Nexp = wf.dot(sum(fc.N)*wf).real
#Nexp = wf.dot(Ntot2*wf).real

print("Density",Nexp)
