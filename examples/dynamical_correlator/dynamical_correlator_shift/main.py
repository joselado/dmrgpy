# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator <Sz_1 delta(w - (H-E0)) Sz_1> of the S=1/2 Heisenberg
# chain, computed with the kernel polynomial method (KPM). The function
# get_distribution accepts a rescaling (a,b) of the energy window. The
# spectrum should not depend on it, which we check for several values of a.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM vectors
sc.nsweeps = 10 # number of DMRG sweeps
sc.kpm_extrapolate = True
sc.kpm_extrapolate_factor = 2.0
e0 = sc.gs_energy()
es = np.linspace(-0.5, 6, 500) # energies
A = sc.Sz[1]

import matplotlib.pyplot as plt
for a in np.linspace(-0.9, 0.5, 4): # shift of the energy window
    (x, y) = sc.get_distribution(X=h-e0, xs=es, delta=1e-2, A=A, B=A, a=a, b=0.9)
    plt.plot(x, y.real, label="a = %.2f" % a)
plt.xlim([-0.5, 4.5])
plt.xlabel("energy w")
plt.ylabel("Dynamical correlator")
plt.legend()
plt.show()
