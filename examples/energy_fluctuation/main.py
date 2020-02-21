# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
import matplotlib.pyplot as plt # library to plot the results

from dmrgpy import spinchain
spins = [3 for i in range(10)] # 2*S+1=2 for S=1/2
for maxm in [1,2,5,10,20,30,40]: # loop over bond dimension
  sc = spinchain.Spin_Hamiltonian(spins) # create spin chain object
  sc.set_exchange(lambda i,j: (abs(i-j)==1)*0.5) # first neighbors
  sc.maxm = maxm # set the bond dimension
  e = sc.gs_energy() # get the ground state energy
  de = sc.gs_energy_fluctuation() # fluctuation
  print("Energy",e,"fluctuation",de,"for bond dimension",maxm)

