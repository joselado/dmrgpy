# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np # conventional numpy library
from dmrgpy import spinchain # library dealing with DMRG for spin chains


####################################
### Create the spin chain object ###
####################################
n = 4 # total number of spins
spins = [4 for i in range(n)] # list with the different spins of your system
# the spins are labeled by 2s+1, so that 2 means s=1/2, 3 meand S=1 ....
sc = spinchain.Spin_Chain(spins) # create the spin chain object



##############################
### Create the hamiltonian ###
##############################
def fj(i,j): # function to define the exchange couplings
    if abs(i-j)==1: return 1.0  # first neighbors
    else: return 0.0 # otherwise
sc.set_exchange(fj) # add the exchange couplings

sc.maxm = 20
o = sc.get_rdm(i=0) # compute density matrix
print(o) # print in a file
print(o.trace())











