# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
import matplotlib.pyplot as plt # library to plot the results

from dmrgpy import spinchain
spins = ["S=1" for i in range(10)] # 2*S+1=2 for S=1/2
sc = spinchain.Spin_Chain(spins) # create spin chain object


h = 0 # initialize
for i in range(len(spins)-1): 
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]


for maxm in [1,2,5,10,20,30,40]: # loop over bond dimension
    sc.set_hamiltonian(h)
    sc.maxm = maxm # set the bond dimension
    e = sc.gs_energy() # get the ground state energy
    de = sc.gs_energy_fluctuation() # fluctuation
    print("Energy",e,"fluctuation",de,"for bond dimension",maxm)










