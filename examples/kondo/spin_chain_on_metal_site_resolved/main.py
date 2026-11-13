# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A chain of 3 spins 1/2 (Heisenberg Js) on a normal metal. The metal is
# a Wilson chain of 5 sites, Kondo coupled (JK) to the spin under the tip.
# We compute the STM spectrum dI/dV with the tip on the end spin and on
# the middle spin. In the ground state of 3 spins the end spin is parallel
# to the total spin and the middle spin is antiparallel to it. So the
# metal screens the chain from the end spin, with peaks close to zero
# bias, and not from the middle spin.
import numpy as np
from dmrgpy import ysr

es = np.linspace(-1.0, 1.0, 201) # bias eV
spectra = []
for site in [2, 1]: # the end spin and the middle spin
    fc = ysr.spin_chain_on_metal(3, 5, site=site, Js=0.3, JK=0.4)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    (e, didv) = ysr.get_stm_spectrum(fc, es, delta=0.03) # tip on site
    print("spin", site, " energy =", round(fc.gs_energy(), 5),
          " dI/dV at zero bias =", round(didv[100], 4))
    spectra.append(didv)

import matplotlib.pyplot as plt
plt.plot(es, spectra[0], linewidth=3, label="tip on the end spin")
plt.plot(es, spectra[1], "--", linewidth=2, label="tip on the middle spin")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
