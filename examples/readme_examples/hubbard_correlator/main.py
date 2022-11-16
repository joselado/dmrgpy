# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import fermionchain
n = 20 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
# first neighbor hopping
h = 0
for i in range(n-1):
  h = h + fc.Cdagup[i]*fc.Cup[i+1]
  h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # Make Hermitian
# Hubbard term
for i in range(n):
  h = h + 4.*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
fc.set_hamiltonian(h) # initialize the Hamiltonian
# compute the two correlators
zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)]
dd = [fc.vev(fc.N[i]).real for i in range(n)]
cc = [fc.vev(fc.Cdagup[0]*fc.Cup[i]).real for i in range(n)]
print(n,len(zz))
print(dd)


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

plt.subplot(1,2,1)
plt.plot(range(len(zz)),zz,marker="o",c="blue",label="$\langle S^z_0S^z_n\\rangle$")
plt.ylabel("Non-local spin correlator")
plt.legend()
plt.xlabel("$n$")
plt.subplot(1,2,2)
plt.plot(range(len(cc)),cc,marker="o",c="red",label="$c^\dagger_0 c_n$")
plt.legend()
plt.ylabel("Non-local charge correlator")
plt.xlabel("$n$")
#plt.xlim([-0.2,4])

plt.tight_layout()
plt.show()





