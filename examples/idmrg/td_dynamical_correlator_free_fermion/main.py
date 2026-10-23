# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')
sys.path.append(os.getcwd()+'/../../../tests') # exact free fermion reference

# Real-time correlator <Sz(x,t) Sz(0,0)> of the dimerized XX chain, which is
# a free fermion chain after a Jordan-Wigner transformation. We compute it
# with iDMRG and a window of sites evolved with TDVP, and compare with the
# exact free fermion result, which is cheap even for thousands of sites.
import numpy as np
from dmrgpy.pyitensor import idmrg, idmrg_window
from _free_fermion_reference import FreeFermionXX

J1, J2 = 1.0, 0.2 # strong and weak bonds
hi = [[J1, ["Sx", 0], ["Sx", 1]], [J1, ["Sy", 0], ["Sy", 1]]] # bond inside the cell
he = [[J2, ["Sx", 1], ["Sx", 2]], [J2, ["Sy", 1], ["Sy", 2]]] # bond between cells
nw = 8 # sites of the window

for seed in range(25): # the solver is random, keep the first usable ground state
    np.random.seed(seed)
    gs = idmrg.idmrg_ground_state([2, 2], hi, he, 2, maxm=24, cutoff=1e-12,
                                  maxiter=150, etol=1e-14, niter=100, verbose=False)
    try:
        idmrg_window.build_window(gs, nw)
        break
    except RuntimeError:
        continue

xs = [-1, 0, 1] # distances
ts, xs, S = idmrg_window.dynamical_correlator_td(gs, n_window=nw, opname_A="Sz",
        opname_B="Sz", dt=0.05, nt=6, cutoff=1e-10, maxdim=40, niter=50,
        x_values=xs, connected=True)
ff = FreeFermionXX(2000, J=J1, J2=J2) # exact reference
x0 = 1000 # reference site, far from the edges
ex = np.array([[ff.sz_sz_connected(x0+int(x), x0, t) for x in xs] for t in ts])
print("Max difference iDMRG-exact =", np.max(np.abs(S-ex)))

import matplotlib.pyplot as plt
for ix, x in enumerate(xs):
    plt.plot(ts, ex[:, ix].real, linewidth=4, color="C%d" % ix, label="exact x=%d" % x)
    plt.plot(ts, S[:, ix].real, "--", linewidth=2, color="k")
plt.xlabel("time t")
plt.ylabel("Re <Sz(x,t) Sz(0,0)>")
plt.legend()
plt.show()
