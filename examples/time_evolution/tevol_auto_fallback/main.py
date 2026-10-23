# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Real-time correlator <Cdag_0(0) N_0(t)> of a fermionic chain with
# tevol_method="AUTO". It uses the fast TEBD propagator when H has only
# nearest-neighbor terms, and falls back to TDVP when it has a longer
# range hopping. We add one hopping between sites 0 and 2 and compare
# AUTO with TDVP directly, which should give the same curve.
import numpy as np
from dmrgpy import fermionchain
from dmrgpy import timedependent

n = 6 # number of sites
def correlator(method, longrange):
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python() # TEBD is in the python backend
    h = 0
    for i in range(n-1): h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    for i in range(n): h = h + 0.3*(-1)**i*fc.N[i] # staggered field
    if longrange: h = h + 0.4*(fc.Cdag[0]*fc.C[2] + fc.Cdag[2]*fc.C[0])
    fc.set_hamiltonian(h)
    fc.tevol_method = method
    return timedependent.evolution_ABA(fc, nt=30, dt=0.05, mode="DMRG",
            A=fc.Cdag[0], B=fc.N[0])

(ts, auto) = correlator("AUTO", True) # falls back to TDVP
(ts, tdvp) = correlator("TDVP", True)
print("Max difference AUTO-TDVP =", np.max(np.abs(np.array(auto)-tdvp)))

import matplotlib.pyplot as plt
plt.plot(ts, np.real(tdvp), linewidth=4, label="TDVP")
plt.plot(ts, np.real(auto), "--", linewidth=2, label="AUTO")
plt.xlabel("time")
plt.ylabel("Re <Cdag_0(0) N_0(t)>")
plt.legend()
plt.show()
