# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Green's function C(t) = <GS|Cup_0(t) Cdagup_0(0)|GS> of the Hubbard
# chain, from real-time evolution. Its Fourier transform is the spectral
# function that ARPES or tunneling would measure. We compare the TEBD
# and TDVP propagators, which should give the same C(t) for short times.
import numpy as np
from dmrgpy import fermionchain
from dmrgpy import timedependent

n = 4 # number of Hubbard sites
t = 1.0 # hopping
U = 4.0 # Hubbard interaction
fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
h = 0
for i in range(n-1):
    h = h + t*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
h = h + h.get_dagger()
for i in range(n):
    h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
fc.setup_python() # TEBD is in the python backend
fc.maxm = 20 # bond dimension
fc.nsweeps = 6 # number of DMRG sweeps
fc.set_hamiltonian(h)

name = (fc.Cup[0], fc.Cdagup[0])
fc.tevol_method = "TEBD"
(ts, c1) = timedependent.evolution_DC(fc, mode="DMRG", name=name, nt=10, dt=0.05)
fc.tevol_method = "TDVP"
(ts, c2) = timedependent.evolution_DC(fc, mode="DMRG", name=name, nt=10, dt=0.05)
print("Max difference TEBD-TDVP =", np.max(np.abs(np.array(c1)-c2)))

import matplotlib.pyplot as plt
plt.plot(ts, np.real(c1), linewidth=4, label="TEBD")
plt.plot(ts, np.real(c2), "--", linewidth=2, label="TDVP")
plt.xlabel("time")
plt.ylabel("Re C(t)")
plt.legend()
plt.show()
