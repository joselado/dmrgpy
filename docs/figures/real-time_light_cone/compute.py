"""Figure 16, real-time_light_cone: a domain wall of the spin-1/2 XX chain
released at t=0, evolved with TDVP, against the exact free-fermion result.

The Jordan-Wigner map of sum_i (Sx_i Sx_i+1 + Sy_i Sy_i+1) is a hopping
1/2, so the dispersion is cos k and the maximal velocity is v=1. The exact
profile is <Sz_i>(t) = C_ii(t) - 1/2, with C(t) = U C(0) U^dagger,
U = exp(-i h t), h the n x n hopping matrix and C(0) diagonal with ones on
the left half.

Run with  bash docs/figures/run.sh real-time_light_cone

Arrays in data.npz:
    ts        (k+1,)    times, t=0 included, in units of 1/J
    sz        (k+1,n)   TDVP <Sz_i>(t), measured on the evolved state
    sz_exact  (k+1,n)   free-fermion <Sz_i>(t)
    x0        ()        position of the wall, between sites n/2-1 and n/2
    params              JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np
import scipy.linalg

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, timedependent

params = common.start(dict(
    model="spin-1/2 XX chain, open, J=1",
    n=40,
    initial_state="ground state of -sum_i eps_i Sz_i, eps=+1 left half, -1 right half",
    tevol_method="TDVP",
    dt=0.05,
    nt_per_chunk=5,
    nchunks=64,
    maxm=100,
    itensor_version=3,
))
n = params["n"]
dt, nt, nch = params["dt"], params["nt_per_chunk"], params["nchunks"]

sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=params["itensor_version"])
sc.maxm = params["maxm"]
sc.tevol_method = params["tevol_method"]
eps = np.array([1.0 if i < n//2 else -1.0 for i in range(n)])
sc.set_hamiltonian(sum(-eps[i]*sc.Sz[i] for i in range(n)))
wf = sc.get_gs()


def profile(wf):
    return np.array([sc.aMb(wf, sc.Sz[i], wf).real for i in range(n)])


sz = [profile(wf)]
print("t=0 profile:", np.round(sz[0], 8), flush=True)
sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                       for i in range(n-1)))

hop = np.zeros((n, n))
for i in range(n-1):
    hop[i, i+1] = hop[i+1, i] = 0.5
c0 = np.diag((eps > 0).astype(float))


def exact(t):
    u = scipy.linalg.expm(-1j*hop*t)
    return np.real(np.diag(u @ c0 @ u.conj().T)) - 0.5


ts = [0.0]
t0 = time.time()
for c in range(1, nch + 1):
    _, _, wf = timedependent.evolve_and_measure(sc, operator=sc.Sz[n//2], nt=nt,
                                                dt=dt, wf=wf, return_wf=True)
    ts.append(c*nt*dt)
    sz.append(profile(wf))
    print("t=%5.2f  max |TDVP-exact| = %.2e  (%.1f s)"
          % (ts[-1], np.abs(sz[-1] - exact(ts[-1])).max(), time.time() - t0),
          flush=True)

ts = np.array(ts)
sz = np.array(sz)
sz_exact = np.array([exact(t) for t in ts])
print("max |TDVP-exact| over all times and sites: %.2e" % np.abs(sz - sz_exact).max())
common.save(__file__, params, ts=ts, sz=sz, sz_exact=sz_exact,
            x0=np.array(n/2 - 0.5))
