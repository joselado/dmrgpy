# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results
####################################
### Create the spin chain object ###
####################################
n = 10 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Chain(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0  # first neighbors
    else: return 0.0 # otherwise

h = 0

for i in range(n):
    for j in range(n): 
        h = h + sc.Sx[i]*sc.Sx[j]*fj(i,j)
        h = h + sc.Sy[i]*sc.Sy[j]*fj(i,j)
        h = h + sc.Sz[i]*sc.Sz[j]*fj(i,j)

for i in range(n):
    h = h + sc.Sx[i]*np.random.random()
    h = h + sc.Sy[i]*np.random.random()
    h = h + sc.Sz[i]*np.random.random()

sc.set_hamiltonian(h)

#sc.itensor_version = "julia"
sc.nsweeps = 1 # one sweep
es = []
wf = None # initial guess
for i in range(10):
    e = sc.gs_energy(wf0=wf,reconverge=True)
    wf = sc.wf0.copy() # save the wavefunction
    es.append(e)
import matplotlib.pyplot as plt
plt.plot(range(len(es)),es,marker="o")
plt.show()











