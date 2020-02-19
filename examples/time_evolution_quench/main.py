# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 2
spins = [2 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain

# create two Hamiltonians
h0 = 0
h1 = 0
h0 = sc.Sz[0] - sc.Sz[1] # AF Hamiltonian
h1 = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1] # Heisenberg


sc.set_hamiltonian(h0) # set AF hamiltonian
sc.get_gs() # compute ground state (it will be stored in the object)
sc.set_hamiltonian(h1) # set the (new) Heisenberg Hamiltonian


from dmrgpy import timedependent
op = sc.Sz[0] # operator to compute
(ts,sz) = timedependent.evolve_and_measure(sc,operator=op,nt=1e4,dt=1e-2)

# now plot the result
import matplotlib.pyplot as plt
plt.plot(ts,sz.real)
plt.show()


