# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import simplechains
# This example shows the quantum phase transition in the trnasverse Ising model
# generate a 1D Ising chain
def get_gap(bx):
    """
    Compute the gap of the 1D Ising model with DMRG
    """
    print("Computing B_x = ",bx)
    b = [bx,0.,0.] # magnetic field
    J = [0.,0.,1.] # Ising coupling
    ssc = simplechains.SSC(s=.5,n=30,b=b,J=J) # generate the Ising chain
    es = ssc.get_excited(mode="DMRG",n=2) # compute the first two energies
    return es[1] - es[0] # return the gap
bs = np.linspace(0.,1.,30) # list of fields
gs = [get_gap(b) for b in bs] # list of gaps
import matplotlib.pyplot as plt
plt.plot(bs,gs,marker="o")
plt.xlabel("Bx")
plt.ylabel("Gap")
plt.show()


