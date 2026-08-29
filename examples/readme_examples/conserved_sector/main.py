# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

from dmrgpy import fermionchain
n = 8 # number of spinless fermionic sites
# quantum-number sectors need itensor_version=3 or "python"
fc = fermionchain.Fermionic_Chain(n,itensor_version="python")
h = 0 # initialize Hamiltonian
for i in range(n-1): h = h + fc.Cdag[i]*fc.C[i+1] # first neighbor hopping
h = h + h.get_dagger() # make the Hamiltonian Hermitian
fc.set_hamiltonian(h) # set the Hamiltonian
nfs = range(n+1) # particle-number sectors
es = [] # storage for the energies
for nf in nfs: # loop over particle-number sectors
  fc.set_conserved_sector(Nf=nf) # confine the calculation to Nf particles
  e = fc.gs_energy() # ground state energy within that sector
  es.append(e) # store
  print("Energy with",nf,"particles",e)
fc.set_conserved_sector() # back to the full Hilbert space
print("Global ground state",fc.gs_energy())


import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 14})

plt.plot(list(nfs),es,marker="o",c="blue")
plt.ylabel("Energy")
plt.xlabel("Number of fermions $N_f$")

plt.tight_layout()
plt.show()
