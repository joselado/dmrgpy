# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A spin S=1 on a superconductor with gap Delta, with a uniaxial
# anisotropy Dz (S^z)^2. One channel can screen only half of the spin, so
# the screened state is a doublet (Sz=+-1/2) and the free state a triplet,
# which Dz splits into Sz=0 and Sz=+-1. We compute the lowest energy of each
# Sz sector against Dz. The lowest one is the ground state, the others are
# the subgap states seen as YSR peaks in dI/dV.
import numpy as np
from dmrgpy import ysr

J, Lambda, n = 0.4, 4.0, 10 # exchange (continuum), discretization, sites
Delta = 0.025 # superconducting gap, the spin is screened at Dz=0
Dzs = np.linspace(-1.0, 1.0, 7)*Delta # anisotropy
E = {0: [], 1: [], 2: []} # energy of the sector 2*Sz
for Dz in Dzs:
    fc = ysr.kondo_chain(n, J=J*ysr.lambda_factor(Lambda), Delta=Delta,
                         Lambda=Lambda, S="1", Dz=Dz)
    fc.maxm = 32 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    r = ysr.subgap_energies(fc, sectors=[0, 1, 2])
    e0 = min(r["sectors"].values()) # ground state energy
    for s in E: E[s].append((r["sectors"][s] - e0)/Delta)
    print("Dz/Delta = %.2f" % (Dz/Delta), r["ground_state"], "2Sz =",
          r["Sz_ground"], "(E-E_GS)/Delta =", [round(E[s][-1], 3) for s in E])

import matplotlib.pyplot as plt
plt.plot(Dzs/Delta, E[0], "o-", label="free, Sz=0")
plt.plot(Dzs/Delta, E[2], "s-", label="free, Sz=+-1")
plt.plot(Dzs/Delta, E[1], "^-", label="screened, Sz=+-1/2")
plt.xlabel("anisotropy Dz/Delta")
plt.ylabel("(E - E_GS)/Delta")
plt.legend()
plt.show()
