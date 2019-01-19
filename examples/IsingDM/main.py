from __future__ import print_function
import sys
import os

import numpy as np
from dmrgpy import simplechains


# This example shows the quantum phase transition in the trnasverse Ising model


# generate a 1D Ising chain
def get_dm(bx):
    """
    Compute the gap of the 1D Ising model with DMRG
    """
    b = [bx,0.,0.] # magnetic field
    J = [0.,0.,1.] # Ising coupling
    ssc = simplechains.SSC(s=.5,n=20,b=b,J=J) # generate the Ising chain
    o = np.sum(ssc.get_dm_dis(n=10)) # return distance between DM
    print("Computing B_x = ",bx,"result",o)
    return o


bs = np.linspace(0.2,2.,5) # list of fields
gs = [get_dm(b) for b in bs] # list of distances

import matplotlib.pyplot as plt
plt.plot(bs,gs,marker="o")
plt.xlabel("Bx")
plt.ylabel("Gap")
plt.show()

