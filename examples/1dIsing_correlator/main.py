# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
# This example shows the quantum phase transition in the transverse Ising model
# generate a 1D Ising chain
def get_gap(bx):
    """
    Compute the gap of the 1D Ising model with DMRG
    """
    print("Computing B_x = ",bx)
    N = 20
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(N)]) # create 
    h = 0
    for i in range(N-1): h = h + sc.Sz[i]*sc.Sz[i+1]
    for i in range(N): h = h + bx*sc.Sx[i]
    sc.set_hamiltonian(h)
    # if you wanted to compute the Mz magnetization instead, uncomment this line
    mz = sum(sc.Sx)/N
    wf = sc.get_gs() # get ground state
    z = wf.dot(mz*wf) # mz
    z2 = wf.dot(mz*(mz*wf)) # mz
    return z2 - z*z

bs = np.linspace(0.,1.,30) # list of fields
gs = [get_gap(b) for b in bs] # list of gaps
import matplotlib.pyplot as plt
plt.plot(bs,gs,marker="o")
plt.xlabel("$B_x$")
plt.ylabel("$\\langle M^2_x \\rangle - \\langle M_x \\rangle^2$")
plt.tight_layout()
plt.show()











