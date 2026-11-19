# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Two spins 1/2 coupled by Js, each on its own metal (its own Wilson chain,
# Kondo coupling JK). We compute the STM spectrum dI/dV with the tip on the
# first spin. With Js = 0 the spin flips at zero bias give one peak. With
# Js > 0 the dimer is a singlet, and the tip excites the triplet: dI/dV
# has two steps at eV = +-Js, the signature measured on spin dimers.
import numpy as np
from dmrgpy import ysr

es = np.linspace(-1.0, 1.0, 201) # bias eV
spectra = []
for Js in [0.0, 0.4]:
    fc = ysr.spin_chain_on_metal(2, 5, site=[0, 1], Js=Js, JK=0.4)
    fc.maxm = 64 # bond dimension
    fc.nsweeps = 16 # number of DMRG sweeps
    # the tip on spin 0, the electron enters its own metal
    (e, didv) = ysr.get_stm_spectrum(fc, es, delta=0.03, imp=fc.impurities[0],
                                     expand_bonds=8, expand_rounds=3)
    print("Js =", Js, " dI/dV at zero bias =", round(didv[100], 4))
    spectra.append(didv)

import matplotlib.pyplot as plt
plt.plot(es, spectra[0], linewidth=3, label="Js = 0")
plt.plot(es, spectra[1], "--", linewidth=2, label="Js = 0.4")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
