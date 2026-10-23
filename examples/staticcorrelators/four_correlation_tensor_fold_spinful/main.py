# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "..", "src"))

# Four-point correlator tensor <C^dag_i C^dag_j C_k C_l> of a spinful
# Hubbard chain in its ground state. The method ctmode="fold" evaluates
# each element as a local operator, and "explicit" builds the full MPO of
# each element. Both give the same tensor, and "fold" is much faster.
import time
import numpy as np
from dmrgpy import fermionchain

times = []
ns = [2, 3, 4] # number of sites
for n in ns:
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n, itensor_version="python")
    h = 0
    for i in range(n-1):
        for f in (0, 1): # spin up and spin down
            a, b = 2*i+f, 2*(i+1)+f
            h = h + fc.Cdag[a]*fc.C[b] + fc.Cdag[b]*fc.C[a] # hopping
    for i in range(n):
        h = h + 2.0*(fc.Cdag[2*i]*fc.C[2*i])*(fc.Cdag[2*i+1]*fc.C[2*i+1]) # Hubbard U=2
    h = h + 0.3*(fc.Cdag[0]*fc.C[0]) # breaks particle-hole symmetry
    fc.set_hamiltonian(h)
    fc.maxm = 20 # bond dimension
    wf = fc.get_gs()
    t0 = time.time()
    t1 = wf.get_four_correlation_tensor(ctmode="fold")
    t2 = time.time()
    t3 = wf.get_four_correlation_tensor(ctmode="explicit")
    times.append([t2-t0, time.time()-t2])
    print("sites", n, "max difference", np.abs(t1-t3).max())
    assert np.abs(t1-t3).max() < 1e-10

import matplotlib.pyplot as plt
times = np.array(times)
plt.plot(2*np.array(ns), times[:, 1], "o-", label="explicit")
plt.plot(2*np.array(ns), times[:, 0], "s-", label="fold")
plt.yscale("log")
plt.xlabel("fermionic modes")
plt.ylabel("time (s)")
plt.legend()
plt.show()
