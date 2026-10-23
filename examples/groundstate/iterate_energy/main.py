# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of a Heisenberg chain in a random magnetic field,
# obtained by repeating one single DMRG sweep. Each call starts from the
# wavefunction of the previous one, so the energy decreases and
# converges to the ground state energy.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# H = sum_i S_i . S_{i+1} + sum_i B_i . S_i, with random B_i
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + np.random.random()*(sc.Sx[i] + sc.Sy[i] + sc.Sz[i])
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.nsweeps = 1 # a single sweep per call

es = []
wf = None # initial guess
for i in range(10):
    e = sc.gs_energy(wf0=wf, reconverge=True)
    wf = sc.wf0.copy() # use this state as the next guess
    es.append(e)

import matplotlib.pyplot as plt
plt.plot(range(len(es)), es, "o-")
plt.xlabel("number of sweeps")
plt.ylabel("energy")
plt.show()
