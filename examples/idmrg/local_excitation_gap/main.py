# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Excitation gap of the transverse-field Ising chain, H = -sum Sx Sx - h sum Sz,
# in the infinite chain (iDMRG). local_excitation_gap() uses the second
# lowest state of a block of the converged iDMRG state, and window=w makes
# the block larger by w sites. The estimate must approach the gap of a
# long finite chain, computed with exact diagonalization (ED).
from dmrgpy import infinitechain, spinchain

J, hz = 1.0, 2.0 # Ising coupling and transverse field, in Pauli units
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "idmrg"
ic.set_hamiltonian(-4*J*ic.SxC[0]*ic.SxR[0] - 2*hz*ic.SzC[0])
ic.maxm = 12 # bond dimension
ic.maxiter = 200
ic.etol = 1e-10
ic.gs_energy()

windows = [0, 1, 2, 3]
gaps = [ic.local_excitation_gap(window=w) for w in windows]

# finite chain of 12 sites as a reference
n = 12
sc = spinchain.Spin_Chain(["1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h - 4*J*sc.Sx[i]*sc.Sx[i+1]
for i in range(n):
    h = h - 2*hz*sc.Sz[i]
sc.set_hamiltonian(h)
gap_ed = sc.get_gap(mode="ED")
print("gaps =", gaps, " ED gap =", gap_ed)

import matplotlib.pyplot as plt
plt.plot(windows, gaps, "o-", label="iDMRG")
plt.axhline(gap_ed, color="gray", linestyle="--", label="ED, 12 sites")
plt.xlabel("window w")
plt.ylabel("excitation gap")
plt.legend()
plt.show()
