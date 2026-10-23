# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spectral function of a non-Hermitian fermion chain,
# a hopping chain with an imaginary onsite potential +i*eta, -i*eta, -i*eta, +i*eta
# repeated every four sites. The correlator <C_i Cdag_i> is computed with
# the biorthogonal left and right eigenstates of H. It gives the
# gain and loss resolved spectrum, as an STM would measure on each site.
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of sites, two cells of four
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1]
h = h + h.get_dagger() # Hermitian hopping
eta = 0.5 # strength of gain and loss
signs = [1, -1, -1, 1] # pattern inside a cell
for i in range(n):
    h = h + 1j*eta*signs[i % 4]*fc.Cdag[i]*fc.C[i]
h = h + 0.03*fc.N[0] # small field, removes a degeneracy
fc.set_hamiltonian(h)

es = np.linspace(-0.1, 3.0, 200) # energies
import matplotlib.pyplot as plt
for i in range(n):
    name = (fc.C[i], fc.Cdag[i])
    (e, d) = fc.get_dynamical_correlator(name=name, mode="ED", es=es, delta=0.3, E_max=10, n=80)
    plt.plot(e, d.real, label="site "+str(i))
plt.xlabel("energy w")
plt.ylabel("spectral function")
plt.legend()
plt.show()
