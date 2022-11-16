# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import fermionchain
import numpy as np
n = 14 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
# first neighbor hopping
h = 0
for i in range(n-1):
  h = h + fc.Cdagup[i]*fc.Cup[i+1]
  h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # Make Hermitian
# Hubbard term
hU = 0
for i in range(n):
  hU = hU + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)

zzs = [] # storage for correlators
Us = np.linspace(0.,4.,6) # Hubbard Us 
for U in Us:
  fc.set_hamiltonian(h+U*hU) # initialize the Hamiltonian
  zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)]
  zzs.append(zz) # store zz correlator


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.rcParams.update({'font.size': 14})

ic = 0
for zz in zzs: 
    c = ic/len(zzs)
    plt.plot(range(len(zz)),zz,marker="o",c=(c,0.,1-c),label="$U="+str(np.round(Us[ic],1))+"$")
    ic +=1
plt.legend()
plt.ylabel("Non-local spin correlator")
plt.xlabel("$n$")

plt.tight_layout()
plt.show()





