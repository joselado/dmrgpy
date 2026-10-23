# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Correlation entropy of the ground state of a chain of spinless fermions
# with interaction U, in the Nambu basis. It is zero for a Slater
# determinant, is largest near U=-1 where the ground state changes,
# and is zero again for strong attraction, where the state is a product.
# We use exact diagonalization, which is enough for four sites.
import numpy as np
from dmrgpy import fermionchain

def get_entropy(U):
    n = 4 # number of sites
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] # hopping
        h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5) # first-neighbor interaction
    fc.set_hamiltonian(h + h.get_dagger())
    wf = fc.get_gs(mode="ED")
    return wf.get_correlation_entropy(basis="Nambu")

Us = np.linspace(-2., 0., 40) # interaction strengths
Ss = [get_entropy(U) for U in Us]

import matplotlib.pyplot as plt
plt.plot(Us, Ss, marker="o")
plt.xlabel("interaction U")
plt.ylabel("correlation entropy")
plt.show()
