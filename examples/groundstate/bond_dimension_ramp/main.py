# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of a Heisenberg-Hubbard chain of spinful fermions with
# DMRG, with and without the bond dimension ramp (bond_ramp=True). With
# the ramp, the first sweeps use a small bond dimension that grows up to
# maxm. The energy is the same and the early sweeps are cheaper.
import time
import numpy as np
from dmrgpy import fermionchain

def energy(maxm, ramp):
    n = 20 # number of sites
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    # hopping, Hubbard U=2 and Heisenberg exchange J=0.4
    h = 0
    for i in range(n-1):
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 2.0*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    for i in range(n-1):
        h = h + 0.4*(fc.Sx[i]*fc.Sx[i+1] + fc.Sy[i]*fc.Sy[i+1] + fc.Sz[i]*fc.Sz[i+1])
    fc.set_hamiltonian(h)
    fc.maxm = maxm # bond dimension
    fc.nsweeps = 12 # number of DMRG sweeps
    fc.bond_ramp = ramp
    t0 = time.time()
    e = fc.gs_energy()
    return e, time.time() - t0

maxms = [20, 40, 60]
flat = [energy(m, False) for m in maxms]
ramp = [energy(m, True) for m in maxms]
print("Max energy difference =", max(abs(f[0]-r[0]) for f, r in zip(flat, ramp)))

import matplotlib.pyplot as plt
plt.plot(maxms, [f[1] for f in flat], "o-", label="flat")
plt.plot(maxms, [r[1] for r in ramp], "s--", label="ramp")
plt.xlabel("bond dimension maxm")
plt.ylabel("time (s)")
plt.legend()
plt.show()
