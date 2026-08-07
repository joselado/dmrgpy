# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import bosonchain
ns = 2 # 6 sites in total, 2 spins
nb = 4 # and 4 bosons
sites = ["S=1/2" for i in range(ns)] # add spins
sites = sites + ["B" for i in range(nb)] # add bosons
bc = bosonchain.SpinBoson_Chain(sites) # create the bosonic chain
h = 0
t = 1.0 # hopping
U = 1.0 # Hubbard
mu = 0.0 # chemical potential
# add the bosonic part
for i in range(ns,ns+nb-1):# loop over bosons (hopping)
    h = h +t*bc.Adag[i]*bc.A[i+1] # hopping
for i in range(ns,ns+nb):# loop over bosons (interaction and chemical potential)
    h = h + U*(bc.N[i]-1.0)*(bc.N[i]-1.0) # onsite interaction
    h = h + mu*bc.N[i] # chemical potential
# add the spin spin part
h = h + bc.Sx[0]*bc.Sx[1]
h = h + bc.Sy[0]*bc.Sy[1]
h = h + bc.Sz[0]*bc.Sz[1]
# and add the spin-boson interaction
h = h + (bc.Sx[0] + bc.Sx[1])*(bc.A[2] + bc.Adag[2])

h = h + h.get_dagger()
bc.set_hamiltonian(h)
e = bc.gs_energy() # compute the ground state energy
print("Energy",e)
navg,n0,n1,n2,n3 = [],[],[],[],[]
for i in range(ns,ns+nb):
    print()
    print("Site #",i)
    navg.append(bc.vev(bc.N[i]).real)
    n0.append(bc.vev(bc.D0[i]).real)
    n1.append(bc.vev(bc.D1[i]).real)
    n2.append(bc.vev(bc.D2[i]).real)
    n3.append(bc.vev(bc.D3[i]).real)
    print("Average occupation",np.round(navg[-1],2))
    print("No occupation",np.round(n0[-1],2))
    print("Single occupation",np.round(n1[-1],2))
    print("Double occupation",np.round(n2[-1],2))
    print("Triple occupation",np.round(n3[-1],2))

import matplotlib.pyplot as plt
bsites = list(range(ns,ns+nb))
plt.plot(bsites,navg,marker="o",label="<N>")
plt.plot(bsites,n0,marker="s",label="P(n=0)")
plt.plot(bsites,n1,marker="^",label="P(n=1)")
plt.plot(bsites,n2,marker="v",label="P(n=2)")
plt.plot(bsites,n3,marker="d",label="P(n=3)")
plt.xlabel("Bosonic site")
plt.ylabel("Occupation / probability")
plt.legend()
plt.show()









