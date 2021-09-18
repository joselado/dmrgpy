# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import parafermionchain
n = 3
pc = parafermionchain.Parafermionic_Chain(n) # create the chain

h = 0
for i in range(n):
  for j in range(n):
      h = h + pc.Sig[i]*pc.Tau[j]*np.random.random()
      h = h + 1j*pc.Sig[i]*pc.Sigd[j]*np.random.random()
      h = h + pc.Tau[i]*pc.Tau[j]*np.random.random()
      h = h + 1j*pc.Tau[i]*pc.Taud[j]*np.random.random()

h = h + h.get_dagger() # Make the Hamiltonian Hermitian
pc.set_hamiltonian(h) # set the Hamiltonian

pc.test_ED()

pc.maxm = 60
pc.nsweeps = 40

ops = pc.Chi + pc.Psi
no = 10

print("Energy ED",pc.gs_energy(mode="ED"))
print("Energy DMRG",pc.gs_energy(mode="DMRG"))

def rando():
    o = ops[np.random.randint(len(ops))]
    return o #+ o.get_dagger()

print(pc.get_excited(mode="ED",n=10))

cs = [rando() for i in range(no)]

c1 = np.array([pc.vev(o,mode="ED") for o in cs])
c2 = np.array([pc.vev(o,mode="DMRG") for o in cs])


plt.plot(range(no),c2.real,color="red",label="DMRG")
plt.scatter(range(no),c1.real,color="blue",label="ED")
plt.legend()
plt.show()









