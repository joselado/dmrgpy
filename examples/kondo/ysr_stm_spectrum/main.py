# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# STM spectrum dI/dV of a spin 1/2 on a superconductor with gap Delta.
# The tip electron tunnels into the substrate under the spin, directly or
# with a spin flip, X^dag = f^dag (U + sigma.S)/2, and with a normal tip
# dI/dV(V) is the spectral function of X at w = eV. Inside the gap we
# see the YSR peaks at +-E_YSR, with different heights at positive and
# negative bias. We compare DMRG with exact diagonalization (ED).
import numpy as np
from dmrgpy import ysr

J, Delta, U = 0.6, 0.1, 0.3 # exchange, gap, direct tunneling
fc = ysr.kondo_chain(4, J=J, Delta=Delta) # spin + 4 Wilson chain sites
fc.maxm = 60 # bond dimension
fc.nsweeps = 15 # number of DMRG sweeps

es = np.linspace(-0.6, 0.6, 241) # bias eV
r = ysr.subgap_energies(fc, mode="ED")
print("E_odd - E_even =", r["E_ysr"], "ground state:", r["ground_state"])
(e1, d1) = ysr.get_stm_spectrum(fc, es, delta=0.01, U=U, mode="DMRG")
(e2, d2) = ysr.get_stm_spectrum(fc, es, delta=0.01, U=U, mode="ED",
                                submode="CVM")
print("Max difference DMRG-ED:", np.max(np.abs(d1 - d2)))
assert np.max(np.abs(d1 - d2)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(e2, d2, linewidth=4, label="ED")
plt.plot(e1, d1, "--", linewidth=2, label="DMRG")
plt.axvline(abs(r["E_ysr"]), color="gray", linestyle=":")
plt.axvline(-abs(r["E_ysr"]), color="gray", linestyle=":")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
