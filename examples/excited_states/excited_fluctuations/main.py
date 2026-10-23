# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Energy fluctuation sqrt(<H^2> - <H>^2) of the lowest excited states of a
# Heisenberg chain with random couplings. For an exact eigenstate it is zero,
# so it measures how well DMRG converged each excited state. We compare
# excited states found with and without Gram-Schmidt orthogonalization.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian with random exchange, H = sum_i J_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + np.random.random()*sc.SS(i,i+1)
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

nex = 3 # number of excited states
def fluctuation(gram_schmidt):
    sc.excited_gram_schmidt = gram_schmidt
    es, wfs = sc.get_excited_states(n=nex)
    hwfs = [h*wf for wf in wfs] # H acting on each state
    return [np.sqrt(abs(hw.dot(hw) - wf.dot(hw)**2)) for hw, wf in zip(hwfs, wfs)]

de0 = fluctuation(True)
de1 = fluctuation(False)

import matplotlib.pyplot as plt
plt.plot(range(nex), de0, "o-", label="with Gram-Schmidt")
plt.plot(range(nex), de1, "o--", label="without Gram-Schmidt")
plt.xlabel("excited state index")
plt.ylabel("energy fluctuation")
plt.legend()
plt.show()
