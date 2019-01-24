# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results
####################################
### Create the spin chain object ###
####################################
n = 40 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0  # first neighbors
    else: return 0.0 # otherwise
sc.set_exchange(fj) # add the exchange couplings
sc.set_fields(lambda i: 3.*np.random.random(3)) # optionally you could add local magnetic fields
sc.nsweeps = 1 # one sweep
es = []
wf0 = None
for i in range(10):
    e = sc.gs_energy(wf0=sc.wf0)
    es.append(e)
import matplotlib.pyplot as plt
plt.plot(range(len(es)),es,marker="o")
plt.show()


