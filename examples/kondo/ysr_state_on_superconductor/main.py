# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A spin 1/2 on a superconductor, coupled by J S.s to the substrate.
# The substrate is a BCS band with gap Delta, mapped to a Wilson chain
# with the pairing Delta on every site. The spin binds a Yu-Shiba-Rusinov
# (YSR) state inside the gap. When T_K is large compared with Delta the
# spin is screened (singlet), when it is small the spin is free (doublet).
# We compute E_odd - E_even, the energy of the subgap peak in dI/dV.
import numpy as np
from dmrgpy import ysr

J, Lambda, n = 0.4, 4.0, 14 # exchange (continuum), discretization, sites
Deltas = np.logspace(-3, -1, 9) # superconducting gaps
es = []
for Delta in Deltas:
    fc = ysr.kondo_chain(n, J=J*ysr.lambda_factor(Lambda), Delta=Delta,
                         Lambda=Lambda)
    fc.maxm = 64 # bond dimension
    fc.nsweeps = 20 # number of DMRG sweeps
    r = ysr.subgap_energies(fc) # ground states of the two parity sectors
    es.append(r["E_ysr"]/Delta)
    print("Delta =", Delta, "(E_odd-E_even)/Delta =", es[-1], r["ground_state"])

import matplotlib.pyplot as plt
plt.semilogx(Deltas, es, "o-")
plt.axhline(0, color="gray", linewidth=0.5)
plt.xlabel("superconducting gap Delta")
plt.ylabel("(E_odd - E_even)/Delta")
plt.title("positive: screened singlet, negative: free doublet")
plt.show()
