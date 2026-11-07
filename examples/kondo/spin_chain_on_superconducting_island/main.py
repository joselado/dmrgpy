# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A chain of 4 spins 1/2 (Heisenberg Js) whose last spin is Kondo coupled
# (JK) to a short superconducting island, 2 sites with the pairing Delta
# put by hand. The island leaks into a normal metal, a Wilson chain of 8
# sites, through the hopping t_c. We compute the STM spectrum dI/dV with
# the tip on the end spin, for a weak and a strong leak into the metal.
# The metal has no gap, so the peaks inside the island gap are resonances,
# and they move to lower bias as t_c grows.
import numpy as np
from dmrgpy import ysr
from dmrgpy.ysrtk import paritystate

es = np.linspace(-1.0, 1.0, 201) # bias eV
spectra = []
for tc in [0.1, 0.4]:
    fc = ysr.island_chain(4, 2, 8, Js=0.3, JK=0.6, Delta=0.25,
                          t_island=0.5, t_c=tc, itensor_version=3)
    fc.maxm = 60 # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    gs, info = paritystate.parity_ground_state(fc) # definite parity
    # nkry=100, the default: lines converged to about 2 percent here,
    # raise nkry for sharp lines
    (e, didv) = ysr.get_stm_spectrum(gs, es, delta=0.03)
    print("t_c =", tc, " island gap =", round(fc.island_gap, 4),
          " E_odd - E_even =", round(info["E_parity"], 3),
          " ground state 2Sz =", info["Sz_ground"])
    spectra.append(didv)

import matplotlib.pyplot as plt
plt.plot(es, spectra[0], linewidth=3, label="t_c = 0.1")
plt.plot(es, spectra[1], "--", linewidth=2, label="t_c = 0.4")
plt.axvline(fc.island_gap, color="gray", linestyle=":")
plt.axvline(-fc.island_gap, color="gray", linestyle=":")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
