# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local dynamical correlator <Sz(t) Sz(0)> of the infinite Heisenberg chain.
# The infinite chain is solved with iDMRG. Then kpm_finite computes the
# correlator in a finite window of n_window cells, with the operators in the
# middle. The result is only reliable if it converges with the window size.
import numpy as np
from dmrgpy import infinitechain

ic = infinitechain.Infinite_Spin_Chain(["1/2"])
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h)
ic.maxm, ic.maxiter, ic.etol = 30, 60, 1e-9 # bond dimension, iterations, tolerance
print("Energy density", ic.gs_energy())

es = np.linspace(-0.5, 4.5, 120) # energies
import matplotlib.pyplot as plt
for n_window in [8, 12, 16]: # size of the window
    (x, y) = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=n_window,
            window_chain_kwargs=dict(maxm=24, nsweeps=10), delta=0.35, es=es)
    plt.plot(es, np.abs(y), label="n_window = %d" % n_window)
plt.xlabel("energy w")
plt.ylabel("|S(w)|")
plt.legend()
plt.show()
