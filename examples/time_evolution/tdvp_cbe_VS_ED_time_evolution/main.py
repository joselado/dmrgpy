# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Quench from the Neel state |up,dn,up,...> to the Heisenberg chain, and
# the magnetization <Sz>(t) of a central site. The Neel state is a product
# state, bond dimension 1. One-site TDVP needs a bond expansion to leave
# it: TDVP_GSE expands only in the first steps, TDVP_CBE (controlled bond
# expansion) at every step. We compare both with exact diagonalization.
import numpy as np
from dmrgpy import spinchain, timedependent

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)], itensor_version="python")
sc.maxm = 32 # bond dimension
sc.nsweeps = 6 # number of sweeps

# the Neel state, as the ground state of a staggered field
sc.set_hamiltonian(sum((-1)**i*sc.Sz[i] for i in range(n)))
wf = sc.get_gs()
wfed = sc.get_gs(mode="ED")

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

op = sc.Sz[n//2] # central spin
z = {}
for method in ["TDVP_GSE", "TDVP_CBE"]:
    sc.tevol_method = method
    ts, z[method] = timedependent.evolve_and_measure(sc, operator=op, wf=wf, nt=60, dt=0.1)
ts, zed = timedependent.evolve_and_measure(sc, operator=op, wf=wfed, nt=60, dt=0.1, mode="ED")
for method in z:
    print("Max difference", method, "- ED =", np.max(np.abs(z[method] - zed)))
assert np.max(np.abs(z["TDVP_CBE"] - zed)) < 1e-4

import matplotlib.pyplot as plt
plt.plot(ts, zed.real, linewidth=4, label="ED")
plt.plot(ts, z["TDVP_GSE"].real, ":", linewidth=2, label="TDVP_GSE")
plt.plot(ts, z["TDVP_CBE"].real, "--", linewidth=2, label="TDVP_CBE")
plt.xlabel("time t")
plt.ylabel("<Sz> of the central site")
plt.legend()
plt.show()
