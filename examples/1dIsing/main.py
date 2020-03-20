# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
# This example shows the quantum phase transition in the trnasverse Ising model
# generate a 1D Ising chain
def get_gap(bx):
    """
    Compute the gap of the 1D Ising model with DMRG
    """
    print("Computing B_x = ",bx)
    sc = spinchain.Spin_Chain([2 for i in range(30)]) # create 
    h = 0
    for i in range(sc.ns-1): h = h + sc.Sz[i]*sc.Sz[i+1]
    for i in range(sc.ns): h = h + bx*sc.Sx[i]
    sc.set_hamiltonian(h)
    return abs(sc.vev(sc.Sz[0]))
    es = sc.get_excited(mode="DMRG",n=2) # compute the first two energies
    return es[1] - es[0] # return the gap
bs = np.linspace(0.,1.,30) # list of fields
gs = [get_gap(b) for b in bs] # list of gaps
import matplotlib.pyplot as plt
plt.plot(bs,gs,marker="o")
plt.xlabel("Bx")
plt.ylabel("Gap")
plt.show()


