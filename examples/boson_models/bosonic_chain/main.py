# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import bosonchain
n = 4

def get_chain(U):
    bc = bosonchain.Bosonic_Chain(n) # create the bosonic chain
    h = 0
    t = 1.0 # hopping
    mu = 0.0 # chemical potential
    for i in range(n-1):  h = h +t*bc.Adag[i]*bc.A[i+1] # hopping
    for i in range(n):  h = h + U*(bc.N[i]-1.0)*(bc.N[i]-1.0) # onsite interaction
    for i in range(n):  h = h + mu*bc.N[i] # chemical potential
    h = h + h.get_dagger()
    bc.set_hamiltonian(h)
    return bc

U0 = 1.0 # Hubbard interaction used for the detailed site-resolved printout
bc = get_chain(U0)
e = bc.gs_energy() # compute the ground state energy
print("Energy",e)
occ_N,occ_D0,occ_D1,occ_D2,occ_D3 = [],[],[],[],[]
for i in range(n):
    print()
    print("Site #",i)
    nv = np.round(bc.vev(bc.N[i]).real,2)
    d0 = np.round(bc.vev(bc.D0[i]).real,2)
    d1 = np.round(bc.vev(bc.D1[i]).real,2)
    d2 = np.round(bc.vev(bc.D2[i]).real,2)
    d3 = np.round(bc.vev(bc.D3[i]).real,2)
    print("Average occupation",nv)
    print("No occupation",d0)
    print("Single occupation",d1)
    print("Double occupation",d2)
    print("Triple occupation",d3)
    occ_N.append(nv); occ_D0.append(d0); occ_D1.append(d1); occ_D2.append(d2); occ_D3.append(d3)
#print("Energy",e)

# also sweep the Hubbard interaction strength U and track the ground
# state energy, a cheap parameter already present in get_chain()
Us = [0.0,0.5,1.0,1.5,2.0,3.0]
Es = []
for U in Us:
    bcU = get_chain(U)
    EU = bcU.gs_energy()
    print("U =",U,"  Energy =",EU)
    Es.append(EU)

import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(10,4))

sites = list(range(n))
axes[0].plot(sites,occ_N,marker="o",label="<N>")
axes[0].plot(sites,occ_D0,marker="s",label="P(n=0)")
axes[0].plot(sites,occ_D1,marker="^",label="P(n=1)")
axes[0].plot(sites,occ_D2,marker="v",label="P(n=2)")
axes[0].plot(sites,occ_D3,marker="d",label="P(n=3)")
axes[0].set_xlabel("Site")
axes[0].set_ylabel("Occupation / probability (U=%s)"%U0)
axes[0].legend()

axes[1].plot(Us,Es,marker="o")
axes[1].set_xlabel("Hubbard U")
axes[1].set_ylabel("Ground state energy")

plt.tight_layout()
plt.show()
