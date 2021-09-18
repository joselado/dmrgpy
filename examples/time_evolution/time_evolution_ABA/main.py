# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 3 # number of sites in the chain
spins = [2 for i in range(n)] # S=1/2
sc = spinchain.Spin_Chain(spins) # create the chain

##################################################################
# This implements shows how to perform the following measurement #
##################################################################
# <X|A^\dagger e^{-iHt} B e^{iHt} A |X> for some state |X>

# create two Hamiltonians
h0 = 0
h1 = 0
for i in range(n): # stagger Hamiltonian
    h0 = h0+(-1)**i*sc.Sz[i]
for i in range(n-1): # Heisenberg Hamiltonian
    h1 = h1+sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1] 

# first get a Neel state by minimizing a dummy Hamiltonian
sc.set_hamiltonian(h0) # set AF hamiltonian
wf = sc.get_gs() # compute ground state
wfED = sc.get_gs(mode="ED") # compute ground state with ED

# now set the new Hamiltonian for the time evolution
sc.set_hamiltonian(h1) # set the (new) Heisenberg Hamiltonian

# define the two operators to apply
A = sc.Sx[0] + sc.Sy[1] # apply some operator to your state
B = sc.Sz[0] # and measure some other operator

from dmrgpy import timedependent

nt = 1e3 # number of time steps
dt = 1e-2 # time step

# perform the two calculations with DMRG and ED
(ts,sz) = timedependent.evolution_ABA(sc,A=A,B=B,nt=nt,dt=dt,wf=wf)
(tsED,szED) = timedependent.evolution_ABA(sc,A=A,B=B,
        nt=nt,dt=dt,wf=wfED,mode="ED")

# now plot the result
import matplotlib.pyplot as plt
plt.plot(ts,sz.real,label="DMRG",c="blue")
plt.scatter(tsED,szED.real,label="ED",c="red")
plt.legend()
plt.show()











