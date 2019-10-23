# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
#sc.set_exchange(lambda i,j: np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])*(abs(i-j)%4==1))
e = sc.gs_energy() # compute the ground state energy
print(e)
# now get the low energy Hamiltonian
h,b = sc.get_effective_hamiltonian()
import scipy.linalg as lg
print(np.round(lg.eigvalsh(b),3))
es,vs = lg.eigh(-b)
vs = vs.transpose()
print(vs[0].real)
print(np.round(lg.eigvalsh(h),3))



