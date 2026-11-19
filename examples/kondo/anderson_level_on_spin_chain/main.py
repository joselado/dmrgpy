# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# A chain of 3 spins 1/2 whose last spin is an Anderson level: an orbital
# with energy ed and repulsion U that hops (tD) onto a metal, a Wilson
# chain of 3 sites. With the tip on the level, dI/dV is its spectral
# function. It has the two Hubbard peaks, at about ed and ed + U, and the
# spin excitations of the chain at low bias. Away from ed = -U/2 the
# spectrum is not symmetric in the bias, since adding and removing an
# electron cost different energies. We compare DMRG with exact
# diagonalization (ED).
import numpy as np
from dmrgpy import ysr
import matplotlib.pyplot as plt

es = np.linspace(-1.5, 1.5, 151) # bias eV
U = 1.6 # repulsion on the level
for ed in [-U/2, -U/4]:
    fc = ysr.anderson_spin_chain(3, 3, Js=0.3, ed=ed, U=U, tD=0.3)
    fc.maxm = 64 # bond dimension
    fc.nsweeps = 16 # number of DMRG sweeps
    (e, d1) = ysr.get_stm_spectrum(fc, es, delta=0.04)
    (e, d2) = ysr.get_stm_spectrum(fc, es, delta=0.04, mode="ED",
                                   submode="CVM")
    print("ed =", ed, " max difference DMRG-ED =", np.max(np.abs(d1 - d2)))
    assert np.max(np.abs(d1 - d2)) < 1e-6
    plt.plot(es, d2, linewidth=4, label="ED, ed = "+str(ed))
    plt.plot(es, d1, "--", linewidth=2, label="DMRG, ed = "+str(ed))
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
