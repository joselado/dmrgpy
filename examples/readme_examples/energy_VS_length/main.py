# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')


from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)] # spins in each site
sc = spinchain.Spin_Chain(spins) # create spin chain object
h = 0 # initialize Hamiltonian
for i in range(len(spins)-1):
  h = h + sc.Sx[i]*sc.Sx[i+1]
  h = h + sc.Sy[i]*sc.Sy[i+1]
  h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h) # create the Hamiltonian


import matplotlib.pyplot as plt





