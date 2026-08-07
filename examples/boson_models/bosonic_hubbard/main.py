# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import bosonchain
n = 3

def get_energy_and_density(U):
    bc = bosonchain.Bosonic_Chain(n) # create the chain

    h = 0 # initialize
    for i in range(n-1): # hopping
        h = h + bc.get_operator("Adag",i)*bc.get_operator("A",i+1)
        h = h + bc.get_operator("Adag",i+1)*bc.get_operator("A",i)

    for i in range(n): # hopping
        den = bc.get_operator("Adag",i)*bc.get_operator("A",i)
        h = h + U*den*den

    bc.set_hamiltonian(h) # setup hamiltonian
    e0 = bc.gs_energy(mode="ED")
    ds = [bc.vev(bc.get_operator("N",i),mode="ED").real for i in range(n)]
    return e0,ds

U0 = 0.2 # interaction strength used in the original example
e0,ds = get_energy_and_density(U0)
print(e0)
print(ds)

# sweep the onsite interaction strength U, a cheap parameter already
# present in this example, and track both the energy and the density
# profile across the chain
Us = [0.0,0.1,0.2,0.4,0.8,1.5]
Es,profiles = [],[]
for U in Us:
    eU,dsU = get_energy_and_density(U)
    print("U =",U,"  Energy =",eU,"  densities =",dsU)
    Es.append(eU)
    profiles.append(dsU)

import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(10,4))

axes[0].plot(range(n),ds,marker="o")
axes[0].set_xlabel("Site")
axes[0].set_ylabel("<N> (U=%s)"%U0)

for U,profile in zip(Us,profiles):
    axes[1].plot(range(n),profile,marker="o",label="U=%s"%U)
axes[1].set_xlabel("Site")
axes[1].set_ylabel("<N>")
axes[1].legend()

plt.tight_layout()
plt.show()
