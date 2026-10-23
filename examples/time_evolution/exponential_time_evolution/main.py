# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Return amplitude <psi|exp(iHt)|psi> of a Heisenberg chain after a spin
# flip on the first site, psi = Sx_0|GS>. It tells how fast the initial
# state is lost under the time evolution of the Hamiltonian. We compute it
# with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.timeevolution import imaginary_exponential # computes exp(iHt)|psi>

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

ts = np.linspace(0.0, 10.0, 50) # times
def overlap(mode):
    wf = sc.get_gs(mode=mode) # ground state
    wf = (sc.Sx[0]*wf).normalize() # flip a spin, no longer an eigenstate
    wfs = imaginary_exponential(h, wf, ts=ts) # exp(iHt)|psi> for every t
    return np.array([wf.dot(w) for w in wfs]).real

u1 = overlap("DMRG")
u2 = overlap("ED")
print("Max difference DMRG-ED =", np.max(np.abs(u1-u2)))

import matplotlib.pyplot as plt
plt.plot(ts, u2, linewidth=4, label="ED")
plt.plot(ts, u1, "--", linewidth=2, label="DMRG")
plt.xlabel("time t")
plt.ylabel("Re <psi|exp(iHt)|psi>")
plt.legend()
plt.show()
