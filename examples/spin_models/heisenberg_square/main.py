# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Frustrated square plaquette of four S=1/2 spins, with nearest neighbor
# coupling J1=1 around the square and a diagonal coupling J2. The ground
# state changes character as J2 grows. The overlap between the ground
# states of two nearby values of J2 drops where this change happens.
import numpy as np
from dmrgpy import spinchain

sc = spinchain.Spin_Chain(["S=1/2" for i in range(4)])
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

h1 = sc.SS(0,1) + sc.SS(1,2) + sc.SS(2,3) + sc.SS(3,0) # square
h2 = sc.SS(0,2) + sc.SS(1,3) # diagonals

j2s = np.linspace(0., 2., 20)
wfs = []
for j2 in j2s:
    sc.set_hamiltonian(h1 + j2*h2)
    wfs.append(sc.get_gs(mode="DMRG")) # ground state
overlaps = [abs(wfs[i-1].dot(wfs[i])) for i in range(1, len(j2s))]

import matplotlib.pyplot as plt
plt.plot((j2s[1:]+j2s[:-1])/2, overlaps, "o-")
plt.xlabel("J2")
plt.ylabel("overlap of neighboring ground states")
plt.show()
