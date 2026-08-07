# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain

n = 20 # number of sites in your chain
spins = ["S=1/2" for i in range(n)] # create the sites
sc = spinchain.Spin_Chain(spins) # create the chain

hfe = sum(sc.Sz) # Fully ferromagnetic
sc.set_hamiltonian(-hfe) ; wfup = sc.get_gs() # Fully up wavefunction
sc.set_hamiltonian(hfe) ; wfdn = sc.get_gs() # Fully down wavefunction

print("Overlap between orthogonal wavefunctions",wfup.dot(wfdn).real)

# Time a repeated local operator application (benchmark only)
wfflip = wfdn.copy()
import time
t0 = time.time()
for i in range(n):
  wfflip = (sc.Sx[0]*sc.Sx[2])*wfflip
#wfflip = (sc.Sx[0])*wfflip
t1 = time.time()
print("Time",t1-t0)

# Now flip all the spins in each site, to go from up to down
wfflip = wfdn.copy() # reset, and perform the actual global spin flip
for i in range(n): wfflip = 2.*sc.Sx[i]*wfflip

# Now do some simple algebra with MPS
print("Rotate and overlap",wfup.dot(wfflip).real)
print("Magnetization up",wfup.dot(hfe*wfup).real)
print("Magnetization down",wfdn.dot(hfe*wfdn).real)

wfs = wfup + 3.7*wfdn # sum the two wavefunctions

print("Overlap with the sum",wfdn.overlap(wfs).real)

# sweep the chain size and see how these quantities evolve with n
def sweep(n):
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    hfe = sum(sc.Sz)
    sc.set_hamiltonian(-hfe) ; wfup = sc.get_gs()
    sc.set_hamiltonian(hfe) ; wfdn = sc.get_gs()
    overlap_updn = wfup.dot(wfdn).real
    wfflip = wfdn.copy()
    for i in range(n): wfflip = 2.*sc.Sx[i]*wfflip
    rotate_overlap = wfup.dot(wfflip).real
    wfs = wfup + 3.7*wfdn
    overlap_sum = wfdn.overlap(wfs).real
    return overlap_updn,rotate_overlap,overlap_sum

ns = range(4,21,4)
out = np.array([sweep(ni) for ni in ns])

import matplotlib.pyplot as plt

plt.plot(ns,out[:,0],marker="o",label="Overlap up/down")
plt.plot(ns,out[:,1],marker="o",label="Rotate and overlap")
plt.plot(ns,out[:,2],marker="o",label="Overlap with sum")
plt.xlabel("Chain size n")
plt.ylabel("Overlap")
plt.legend()
plt.show()










