# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# dI/dV of a spin S=1 with anisotropy Dz<0 on a superconductor. The ground
# state is a pair Sz=+-1, and we solve it inside one Sz sector. TASK works
# in the tangent space of this state, and its small bonds make this space
# too small: the spectrum is off. pad_bonds="full" pads every bond with
# empty directions, the state does not change but the tangent space grows
# to the whole Hilbert space, and the spectrum is the exact one.
import numpy as np
from dmrgpy import ysr
from dmrgpy.ysrtk import paritystate

fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, S="1", Dz=-0.05)
fc.maxm = 64 # bond dimension
fc.nsweeps = 15 # number of DMRG sweeps
es = np.linspace(-1.2, 1.2, 241) # bias eV
_, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, U=0.3, mode="ED",
                             submode="CVM") # exact reference
c, info = paritystate.parity_ground_state(fc) # state of one Sz sector
_, plain = ysr.get_stm_spectrum(c, es, delta=0.08, U=0.3) # TASK
_, padded = ysr.get_stm_spectrum(c, es, delta=0.08, U=0.3, pad_bonds="full")
print("Max difference without padding", np.max(np.abs(plain - ed)))
print("Max difference with padding", np.max(np.abs(padded - ed)))
assert np.max(np.abs(padded - ed)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(es, ed, linewidth=4, label="ED")
plt.plot(es, plain, "--", label="TASK")
plt.plot(es, padded, ":", label="TASK, pad_bonds=\"full\"")
plt.xlabel("bias eV")
plt.ylabel("dI/dV")
plt.legend()
plt.show()
