# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import parafermionchain
n = 3
pc = parafermionchain.Parafermionic_Chain(n,Z=4) # create the chain
pc.test_ED()

h = 0
for i in range(n):
  for j in range(n):
      h = h + pc.N[i]*pc.N[j]*np.random.random()
      h = h + pc.Sig[i]*pc.Sig[j]*np.random.random()
      h = h + 1j*pc.Sig[i]*pc.Sigd[j]*np.random.random()
      h = h + pc.Tau[i]*pc.Tau[j]*np.random.random()
      h = h + 1j*pc.Tau[i]*pc.Taud[j]*np.random.random()

h = h + h.get_dagger() # Make the Hamiltonian Hermitian
pc.set_hamiltonian(h) # set the Hamiltonian
pc.maxm = 60
pc.kpmmaxm = 60
print(pc.gs_energy(mode="DMRG"))
print(pc.gs_energy(mode="ED"))

print(pc.get_excited(mode="ED",n=10))

ops = pc.Chi + pc.Psi
no = 10

def rando():
    o = ops[np.random.randint(len(ops))]
    return o #+ o.get_dagger()

name = [rando(),rando()]
delta = 1e-1
(e1,d1) = pc.get_dynamical_correlator(name=name,mode="ED",submode="INV",
        delta=delta)
(e2,d2) = pc.get_dynamical_correlator(name=name,mode="ED",submode="KPM",
        delta=delta)
(e0,d0) = pc.get_dynamical_correlator(name=name,mode="DMRG",submode="KPM",
        delta=delta)

plt.plot(e0,d0.real,color="red",label="DMRG")
plt.scatter(e2,d2.real,color="green",label="ED-KPM")
plt.plot(e1,d1.real,color="blue",label="ED-INV")
plt.legend()
plt.show()
