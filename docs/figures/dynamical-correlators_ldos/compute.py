"""Figure 17, dynamical-correlators_ldos: local density of states A_i(omega) of
an open spinless fermion chain at half filling, at a bulk site and at the
edge site, for V=0 (free fermions, exact reference) and V=3t (charge
density wave, gapped), from get_spectral_function (SECTOR method).

H = -t sum_i (c^dag_i c_{i+1} + h.c.) + V sum_i (n_i - 1/2)(n_{i+1} - 1/2)

Run with  bash docs/figures/run.sh dynamical-correlators_ldos

Arrays in data.npz:
    Vs          (nv,)         interaction strengths V, in t
    sites       (ns,)         sites i of A_i (bulk, edge)
    es          (ne,)         omega - mu, in t
    A           (nv, ns, ne)  A_i(omega - mu) from SECTOR, real part
    captured    (nv, ns, 2)   captured weight, particle and hole side
    mu          (nv,)         chemical potential (E0(N+1)-E0(N-1))/2, in t
    gap         (nv,)         charge gap E0(N+1)+E0(N-1)-2E0(N), in t
    A_free      (ns, ne)      exact free-fermion A_i at V=0, same delta
    eps         (n,)          free-fermion single-particle levels, in t
    gap_free    ()            free-fermion gap eps_{N+1}-eps_N, in t
    params                    JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import fermionchain

params = common.start(dict(
    model="spinless fermions, open chain, hopping -t, t=1, "
          "V sum (n_i-1/2)(n_i+1-1/2), half filling",
    n=int(os.environ.get("FIG_N", 20)),
    Vs=[0.0, 3.0],
    maxm=80,
    nsweeps=12,
    cutoff=1e-12,
    nex=int(os.environ.get("FIG_NEX", 50)),
    delta=0.1,
    es=[-5.0, 5.0, 1000],
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]
N = n // 2
sites = [n // 2, 0]
delta = params["delta"]
es = np.linspace(*params["es"][:2], params["es"][2])


def lorentz(x, x0):
    return (delta/np.pi) / ((x - x0)**2 + delta**2)


# exact free-fermion reference on the same open chain
hop = -(np.eye(n, k=1) + np.eye(n, k=-1))
eps, phi = np.linalg.eigh(hop)
mu_free = 0.5*(eps[N-1] + eps[N])
A_free = np.array([sum(phi[i, k]**2*lorentz(es, eps[k] - mu_free)
                       for k in range(n)) for i in sites])
print("free fermions: mu = %.3e, gap = %.5f" % (mu_free, eps[N] - eps[N-1]))

A, cap, mus, gaps = [], [], [], []
for V in params["Vs"]:
    fc = fermionchain.Fermionic_Chain(n, itensor_version=params["itensor_version"])
    h = 0
    for i in range(n-1):
        h = h - fc.Cdag[i]*fc.C[i+1] - fc.Cdag[i+1]*fc.C[i]
        h = h + V*(fc.N[i] - 0.5)*(fc.N[i+1] - 0.5)
    fc.set_hamiltonian(h)
    fc.maxm = params["maxm"]
    fc.nsweeps = params["nsweeps"]
    fc.cutoff = params["cutoff"]
    fc.set_conserved_sector(Nf=N)
    Av, cv = [], []
    for i in sites:   # the same chain, so the sector states are reused
        t0 = time.time()
        _, a, info = fc.get_spectral_function(i=i, nex=params["nex"], es=es,
                                              delta=delta, return_poles=True)
        Av.append(np.asarray(a).real)
        cv.append([info["particle"]["captured"], info["hole"]["captured"]])
        print("V=%.1f i=%2d  mu=%.5f  gap=%.5f  captured particle %.4f "
              "hole %.4f  (%.1f s)" % (V, i, info["mu"], info["gap"],
                                       cv[-1][0], cv[-1][1], time.time()-t0),
              flush=True)
        if V == 0.0:
            print("        max|A-A_free| / max A_free = %.3e"
                  % (np.max(np.abs(Av[-1] - A_free[len(Av)-1]))
                     / np.max(A_free[len(Av)-1])), flush=True)
    A.append(Av); cap.append(cv); mus.append(info["mu"]); gaps.append(info["gap"])

common.save(__file__, params, Vs=np.array(params["Vs"]), sites=np.array(sites),
            es=es, A=np.array(A), captured=np.array(cap), mu=np.array(mus),
            gap=np.array(gaps), A_free=A_free, eps=eps,
            gap_free=np.array(eps[N] - eps[N-1]))
