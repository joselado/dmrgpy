# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the chain


def get(U):
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
        h = h + U*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    
    h = h + h.get_dagger()
    
    fc.set_hamiltonian(h) # Hamiltonian
    wf = fc.get_gs(mode="ED") # compute ground state
    h = 0
    
    dm = wf.get_dm(inds=range(n))
    import scipy.linalg as lg
    return lg.eigvalsh(dm)


import matplotlib.pyplot as plt

Us = np.linspace(0.,20.0,20)

for U in Us:
    y = get(U)
    plt.scatter(0*y+U,y)

plt.show()










