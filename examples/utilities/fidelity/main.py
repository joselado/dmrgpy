# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fidelity susceptibility of the transverse field Ising chain,
# H = 4 sum_i Sz_i Sz_{i+1} + l 2 sum_i Sx_i.
# It measures how fast the ground state changes when the field l changes,
# and it has a peak at the quantum phase transition.
import numpy as np
from dmrgpy import spinchain, fidelity

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h0 = 0 # Ising coupling, in a ring
for i in range(n):
    h0 = h0 + 4*sc.Sz[i]*sc.Sz[(i+1) % n]
h1 = 0 # transverse field
for i in range(n):
    h1 = h1 + 2*sc.Sx[i]

ls = np.linspace(0., 2., 10) # field strengths
fs = [fidelity.get_fidelity(sc, h0, h1, l, n=3, mode="ED") for l in ls]

import matplotlib.pyplot as plt
plt.plot(ls, fs, "o-")
plt.xlabel("transverse field l")
plt.ylabel("fidelity susceptibility")
plt.show()
