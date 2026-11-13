# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# STM spectrum of a chain of three spins 1/2 on a metal, with the tip on
# the end spin, which is Kondo coupled to the metal. The numerical
# renormalization group (NRG) starts from the bond basis of a DMRG state
# of the spins and the first Wilson sites (the head), and continues the
# Wilson chain from there (the tail). We compare it with the exact NRG,
# whose first shell holds all 8 states of the spins.
import numpy as np
from dmrgpy import ysr

Lambda = 4.0 # discretization of the Wilson chain
JK = 0.3*ysr.lambda_factor(Lambda) # rho JK = 0.15 in the continuum
fc = ysr.spin_chain_on_metal(3, 12, Js=0.01, JK=JK, Lambda=Lambda,
                             itensor_version=3) # 3 spins, 12 Wilson sites
fc.maxm, fc.nsweeps = 160, 30 # bond dimension and sweeps of the head

es = np.logspace(-8, 0, 50) # positive biases
(e, A1) = ysr.get_stm_spectrum(fc, es, submode="NRG", head="exact", keep=1600)
(e, A2) = ysr.get_stm_spectrum(fc, es, submode="NRG", head="dmrg", keep=1600,
                               m=4, nkry=30) # hand over 4 sites before the end
low = es < 1e-3 # below the scale of the hand-over, where the tail is used
print("Max relative difference below 1e-3:", np.max(np.abs(A2/A1-1)[low]))
assert np.max(np.abs(A2/A1-1)[low]) < 0.03

import matplotlib.pyplot as plt
plt.semilogx(es, A1, linewidth=4, label="exact first shell")
plt.semilogx(es, A2, "--", linewidth=2, label="DMRG head, NRG tail")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
