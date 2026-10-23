# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlation entropy after creating an electron-hole pair on top of the
# ground state of spinless fermions. It is zero for a Slater determinant,
# so after a quench it measures how interactions V build up correlations.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
nt = 100 # number of time steps
dt = 0.1 # time step
import matplotlib.pyplot as plt
for V in [0., 1., 2.]: # nearest-neighbor interaction
    fc = fermionchain.Fermionic_Chain(n)
    h0, hv = 0, 0
    for i in range(n-1):
        h0 = h0 + fc.Cdag[i]*fc.C[i+1] # hopping
        hv = hv + (fc.N[i]-.5)*(fc.N[i+1]-.5) # interaction
    h = h0 + h0.get_dagger() + V*hv
    fc.set_hamiltonian(h)
    wf = fc.get_gs(mode="ED")
    wf = (fc.Cdag[1]*fc.C[0]*wf).normalize() # electron-hole pair
    ents = []
    for i in range(nt): # time evolution with e^{iHdt}
        ents.append(wf.get_correlation_entropy())
        wf = fc.exponential(1j*dt*h, wf)
    plt.plot(dt*np.arange(nt), ents, label="V = "+str(V))
plt.xlabel("time")
plt.ylabel("correlation entropy")
plt.legend()
plt.show()
