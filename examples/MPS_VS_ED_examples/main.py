# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(2)] # or just a dimer with S=1/2
sc = spinchain.Spin_Chain(spins) # create the spin chain object

H_I = sc.Sz[0]*sc.Sz[1] # this would be an Ising Hamiltonian
H_XY = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] # this would be an XY Hamiltonian
H_Hei = sc.Sx[0]*sc.Sx[1] + sc.Sy[0]*sc.Sy[1] + sc.Sz[0]*sc.Sz[1] # this would be a Heisenberg Hamiltonian


sc.set_hamiltonian(H_I) # set the Hamiltonian
E0_I = sc.gs_energy(mode="ED") # compute the energy for the Ising model

sc.set_hamiltonian(H_XY) # set the Hamiltonian
E0_XY = sc.gs_energy(mode="ED") # compute the energy for the Ising model

sc.set_hamiltonian(H_Hei) # set the Hamiltonian
E0_Hei = sc.gs_energy(mode="ED") # compute the energy for the Ising model
print("Energy of the Ising model",E0_I)
print("Energy of the XY model",E0_XY)
print("Energy of the Heisenberg model",E0_Hei)
