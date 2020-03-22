# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results
####################################
### Create the spin chain object ###
####################################
n = 8 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Chain(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: 
        return np.array([[1,0,0],[0,1,0],[0,0,1]])  # first neighbors
    else: return 0.0 # otherwise
sc.set_exchange(fj) # add the exchange couplings
#sc.set_fields(lambda i: np.random.random(3)) # optionally you could add local magnetic fields
# parameters controlling the DMRG algorithm, in principle default ones are fine
#sc.maxm = 40 # maximum bond dimension
#sc.nsweeps = 12 # number of DMRG sweeps
############################################################
# Perform ground state calculation and compute correlators #
############################################################
# compute ground state energy
sc.cutoff = 1e-10
#e0 = sc.gs_energy() # compute ground state energy
#print("Energy",e0/sc.ns)
#sc.get_gs(n=7)
# this array constains the pairs of spins on which you want to compute
# <GS|S_i S_j GS>
pairs = [(0,i) for i in range(n)] # between the edge and the rest
cs = sc.get_correlator(pairs=pairs,mode="DMRG").real # get the static correlators
cs1 = sc.get_correlator(pairs=pairs,mode="ED").real # get the static correlators
########################
# Now plot the results #
########################
# now plot the result
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(range(n),cs,marker="o",c="blue",label="DMRG")
plt.scatter(range(n),cs1,marker="o",c="red",s=100,label="ED")
plt.legend()
plt.xlabel("N")
plt.ylabel("<GS|S_0 S_N |GS>")
plt.show()


