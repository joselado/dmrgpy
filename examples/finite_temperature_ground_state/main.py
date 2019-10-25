# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains
import matplotlib.pyplot as plt # library to plot the results
####################################
### Create the spin chain object ###
####################################
n = 2 # total number of spins
spins = [2 for i in range(n)] # list with the different spins of your system
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain object
##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: 
        return np.array([[1,0,0],[0,1,0],[0,0,1]])  # first neighbors
    else: return 0.0 # otherwise
sc.set_exchange(fj) # add the exchange couplings
#sc.set_fields(lambda i: [i==0,0.,0.]) 
############################################################
# Perform ground state calculation and compute correlators #
############################################################
Ts = np.linspace(0.0,1.0,10) # temperatures
def f(T):
    return sc.gs_energy(T=T,mode="ED")
es = [f(T) for T in Ts] # energies

########################
# Now plot the results #
########################
# now plot the result
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(Ts,es,marker="o",c="blue",label="ED")
#plt.scatter(range(n),cs1,marker="o",c="red",s=100,label="ED")
plt.legend()
plt.xlabel("N")
plt.ylabel("Energy")
plt.show()


