# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Imaginary-time evolution e^{-tau H}|Neel> of the S=1/2 Heisenberg chain,
# by implicit trapezoid steps on a logarithmic grid (arXiv:2606.02930).
# The steps double every few steps, so tau=1000 takes 72 steps.
# The energy relaxes to the ground-state energy; we check it against
# exact diagonalization (ED) on the same grid.
import numpy as np
from dmrgpy import spinchain

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

# the Neel state, as the ground state of a staggered field
sc.set_hamiltonian(sum((-1)**i*sc.Sz[i] for i in range(n)))
wf = sc.get_gs()
wfed = sc.get_gs(mode="ED")

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

tau = 1000. # final imaginary time
taus, e = sc.evolve_imaginary_time_and_measure(h, wf, tau, nsteps=4)
taus, eed = sc.evolve_imaginary_time_and_measure(h, wfed, tau, mode="ED")
e0 = sc.gs_energy(mode="ED") # ground-state energy
print("Number of steps =", len(taus)-1)
print("Max difference DMRG-ED =", np.max(np.abs(e-eed)))
assert np.max(np.abs(e-eed)) < 1e-2

import matplotlib.pyplot as plt
plt.plot(taus[1:], (eed-e0).real[1:], linewidth=4, label="ED")
plt.plot(taus[1:], (e-e0).real[1:], "o--", label="DMRG, log grid")
plt.xscale("log") ; plt.yscale("log")
plt.xlabel("imaginary time tau")
plt.ylabel("E(tau) - E0")
plt.legend()
plt.show()
