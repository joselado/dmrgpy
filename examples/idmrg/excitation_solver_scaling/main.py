# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Excitation energies E(k) of the transverse field Ising chain, from the
# tangent-space ansatz on top of a VUMPS ground state. The eigenproblem
# can be solved by building the matrix (dense) or with Lanczos (iterative).
# We time both against the bond dimension D and check that they agree.
# H = -sigma^x_i sigma^x_i+1 - g sigma^z_i, gapped for g > 1.
import time
import numpy as np
from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg_excitations

g = 2.5
ks = np.linspace(0, np.pi, 7) # momenta

def scan(env, threshold): # threshold forces dense (large) or iterative (0)
    idmrg_excitations._DENSE_EIG_MAX = threshold
    env.resolvent_cache.clear()
    t0 = time.time()
    e = [idmrg_excitations.excitation_energies(env, k, n=1)[0].real for k in ks]
    return time.time()-t0, np.array(e)

Ds = [2, 3, 4, 6]
td, ti, diff = [], [], []
for D in Ds:
    ic = infinitechain.Infinite_Spin_Chain(["1/2"])
    ic.gs_method = "vumps"
    ic.maxm = D # bond dimension
    ic.etol = 1e-12
    ic.vumps_nrestarts = 4
    ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
    ic.gs_energy()
    env = ic._get_excitation_environment()
    t1, e1 = scan(env, 10**9)
    t2, e2 = scan(env, 0)
    td.append(t1) ; ti.append(t2) ; diff.append(np.max(np.abs(e1-e2)))
print("Max difference dense-iterative =", max(diff))

import matplotlib.pyplot as plt
plt.plot(Ds, td, "o-", linewidth=4, label="dense")
plt.plot(Ds, ti, "s--", linewidth=2, label="iterative")
plt.xlabel("bond dimension D")
plt.ylabel("time of the k scan (s)")
plt.legend()
plt.show()
