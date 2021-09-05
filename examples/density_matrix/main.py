# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')
#import os ; import sys ; sys.path.append(os.environ["PYGRAROOT"])

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Chain(n) # create the chain

h = 0
U = 0.3

for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + U*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)


h = h + h.get_dagger()

fc.set_hamiltonian(h) # Hamiltonian
wf = fc.get_gs(mode="ED") # compute ground state

dm = wf.get_dm(inds=range(n))
#print(dm)
import scipy.linalg as lg
print(np.round(lg.eigvalsh(dm),3))
