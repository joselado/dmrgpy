# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Correlation entropy of the ground state of an interacting spinless fermion chain,
# S = -sum_k [p_k log p_k + (1-p_k) log(1-p_k)], with p_k the eigenvalues of <Cdag_i C_j>.
# It is zero for a Slater determinant and grows with the interaction U,
# so it measures how far the state is from a non-interacting one.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites

def entropy(U):
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] # hopping
        h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5) # interaction
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    wf = fc.get_gs(mode="ED")
    return wf.get_correlation_entropy(mode="ED")

Us = np.linspace(-2, 0, 30)
Ss = [entropy(U) for U in Us]

import matplotlib.pyplot as plt
plt.plot(Us, Ss)
plt.xlabel("interaction U")
plt.ylabel("correlation entropy")
plt.show()
