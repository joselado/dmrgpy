# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Four-point correlator tensor <Cdag_i C_j Cdag_k C_l> of a fermion chain with
# nearest-neighbor interaction, computed from the ground state with two algorithms.
# ctmode="sweep" contracts one index tuple at a time, ctmode="batched"
# groups the tuples that share partial contractions, so it is faster.
# Both give the same tensor, so we print the maximum difference and the times.
import time
import numpy as np
from dmrgpy import fermionchain

def get_gs(n):
    fc = fermionchain.Fermionic_Chain(n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
        h = h + fc.N[i]*fc.N[i+1] # interaction
    fc.set_hamiltonian(h)
    fc.maxm = 20 # bond dimension
    fc.nsweeps = 4 # number of DMRG sweeps
    return fc.get_gs()

ns = [4, 6, 8, 10] # chain lengths
t_batched, t_sweep, diff = [], [], []
for n in ns:
    wf = get_gs(n)
    t0 = time.time()
    ct1 = wf.get_four_correlation_tensor(ctmode="batched")
    t1 = time.time()
    ct2 = wf.get_four_correlation_tensor(ctmode="sweep")
    t2 = time.time()
    t_batched.append(t1-t0) ; t_sweep.append(t2-t1)
    diff.append(np.max(np.abs(ct1-ct2)))
print("Max difference batched-sweep =", max(diff))
assert max(diff) < 1e-10

import matplotlib.pyplot as plt
plt.plot(ns, t_sweep, "o-", label="sweep")
plt.plot(ns, t_batched, "s--", label="batched")
plt.xlabel("number of sites")
plt.ylabel("time [s]")
plt.legend()
plt.show()
