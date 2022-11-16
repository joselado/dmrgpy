# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')


from dmrgpy import spinchain
n = 20 ; spins = ["S=1" for i in range(n)] # S=1 chain
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
(e0,d0) = sc.get_dynamical_correlator(name=(sc.Sz[0],sc.Sz[0]))
(eb,db) = sc.get_dynamical_correlator(name=(sc.Sz[n//2],sc.Sz[n//2]))


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

plt.plot(e0,d0,marker="o",c="blue",label="edge")
plt.plot(eb,db,marker="o",c="red",label="bulk")
plt.legend()
plt.ylabel("$\langle S^z_n \delta(\omega- H +E_0) S^z_{n} \\rangle$")
plt.xlabel("$\omega$")
plt.xlim([-0.2,4])

plt.tight_layout()
plt.show()





