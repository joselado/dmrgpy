# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Thermal energy E(T) = Tr(H exp(-H/T))/Z of the S=1/2 Heisenberg chain.
# METTS samples product states and evolves each one to beta/2 = 1/(2T),
# here by implicit steps on a logarithmic grid (arXiv:2606.02930), on the
# Julia backend. At low T the energy stays close to the ground-state
# energy, and it rises as T approaches the spin gap. We compare with
# exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.setup_julia() # ITensors.jl backend
sc.maxm = 20 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps, for the ground-state energy

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

Ts = [0.2, 0.3, 0.5, 0.7] # temperatures
ed = [sc.vev(h, mode="ED", T=T).real for T in Ts] # exact thermal energy
metts, err = [], []
for T in Ts:
    (e, de) = sc.metts_vev(h, T, nsamples=50, nwarmup=5, itime="implicit", seed=1)
    metts.append(e.real) ; err.append(de)
print("Max difference METTS-ED =", np.max(np.abs(np.array(metts)-ed)))
# the error bar grows with T (0.02 at T=0.2, 0.1 at T=0.7 with 50 samples),
# so we check every point against four of its own error bars
assert np.all(np.abs(np.array(metts)-ed) < 4*np.array(err))

import matplotlib.pyplot as plt
Tl = np.linspace(0.1, 0.8, 40) # temperatures of the exact curve
plt.plot(Tl, [sc.vev(h, mode="ED", T=T).real for T in Tl], linewidth=4, label="ED")
plt.errorbar(Ts, metts, yerr=err, fmt="o--", label="METTS, implicit steps")
plt.xlabel("temperature T")
plt.ylabel("energy E(T)")
plt.legend()
plt.show()
