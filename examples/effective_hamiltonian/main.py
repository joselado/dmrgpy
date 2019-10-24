# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
e = sc.gs_energy() # compute the ground state energy
# now get the low energy Hamiltonian
h,b = sc.get_effective_hamiltonian()
print(sc.get_excited(n=4))
import scipy.linalg as lg
print(np.round(lg.eigvalsh(b),3))
es,vs = lg.eigh(-b)
vs = vs.transpose()
print(vs[0].real)
print(np.round(lg.eigvalsh(h),3))



