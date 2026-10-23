# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin and charge correlators in the ground state of the Hubbard chain at
# half filling, with U/t=4. We compute <Sz_0 Sz_j> and <c^dagger_0 c_j>.
# The spin correlator is antiferromagnetic, and both decay with distance.
from dmrgpy import fermionchain

n = 20 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger() # make it Hermitian
for i in range(n):
    h = h + 4.*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5) # Hubbard term
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)] # spin correlator
cc = [fc.vev(fc.Cdagup[0]*fc.Cup[i]).real for i in range(n)] # hopping correlator

import matplotlib.pyplot as plt
plt.plot(range(n), zz, "o-", label="<Sz_0 Sz_j>")
plt.plot(range(n), cc, "s-", label="<c+_0 c_j>")
plt.xlabel("site j")
plt.ylabel("correlator")
plt.legend()
plt.show()
