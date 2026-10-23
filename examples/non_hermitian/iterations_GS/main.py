# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of a non-Hermitian Heisenberg chain with an imaginary field,
# H = sum_i S_i . S_{i+1} + 0.4 i sum_i Sz_i.
# The MPS ground state is found with an Arnoldi iteration. We
# see how the error of the energy goes down with the number of iterations.
import numpy as np
from dmrgpy import spinchain

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + 1j*0.4*sc.Sz[i] # non-Hermitian part
sc.set_hamiltonian(h)
sc.maxm = 30 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
e0 = sc.gs_energy(mode="ED") # exact energy

ites = [1, 4, 10, 40] # number of Arnoldi iterations
des = []
for ite in ites:
    sc.set_hamiltonian(h) # start again for each ite
    e = sc.gs_energy(mode="DMRG", maxit=ite)
    des.append(np.abs(e0-e))
    print("iterations =", ite, "error =", des[-1])

import matplotlib.pyplot as plt
plt.plot(ites, des, "o-")
plt.yscale("log")
plt.xlabel("Arnoldi iterations")
plt.ylabel("energy error")
plt.show()
