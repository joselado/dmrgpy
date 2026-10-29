# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Two-level Anderson impurity, each level coupled to its own Wilson chain,
# with the NRG-MPS method of arXiv:2501.10746. The repulsion U acts between
# all four spin-orbitals, so the impurity holds one electron. When the two
# couplings are equal the ground state is an SU(4) Kondo singlet, and the
# four Schmidt coefficients of the singly occupied impurity are all 1/2.
# When they differ the coefficients split in two pairs, and on a longer
# chain (lower energy) they go toward 1/sqrt(2), 1/sqrt(2), 0, 0: the SU(2)
# Kondo singlet of the level with the larger coupling.
import numpy as np
from dmrgpy import wilsonchain

U, t1, Lambda = 1.0, 0.3, 3.0 # repulsion, coupling of level 1, discretization
n = 7 # Wilson chain sites for each level, odd

t2s = np.linspace(0.15, 0.45, 7) # couplings of level 2
lam, M = [], []
for t2 in t2s:
    fc = wilsonchain.anderson_chain(n, ed=[-U/2, -U/2], U=U, tD=[t1, t2],
                                    Lambda=Lambda)
    fc.maxm = 100 # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    fc.cutoff = 1e-12
    r = fc.get_kondo_screening()
    lam.append(r["schmidt"])
    M.append(r["M_imp"]) # single occupancy of each level
    print("t2 =", t2, "Schmidt coefficients", r["schmidt"])

import matplotlib.pyplot as plt
plt.plot(t2s/t1, lam, "o-")
plt.plot(t2s/t1, np.array(M)[:, 0], "k--", label="M_1")
plt.plot(t2s/t1, np.array(M)[:, 1], "k:", label="M_2")
plt.axhline(0.5, color="gray", linewidth=0.5)
plt.xlabel("t2/t1")
plt.ylabel("Schmidt coefficients and M_l")
plt.legend()
plt.show()
