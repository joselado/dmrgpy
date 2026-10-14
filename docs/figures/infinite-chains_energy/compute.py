"""Figure 20, infinite-chains_energy: energy per site of the infinite
spin-1/2 Heisenberg chain against the bond dimension (iDMRG), and the
spin correlator at long distance.

Run with  bash docs/figures/d22.sh run infinite-chains_energy

Arrays in data.npz:
    maxm       (k,)  bond dimensions
    energy     (k,)  iDMRG energy per site e(maxm), in units of J
    converged  (k,)  ic.converged after each solve (bool)
    time_s     (k,)  wall time of each solve, in s
    e_exact    ()    Bethe-ansatz energy per site 1/4 - ln 2
    r          (m,)  distances of the correlator
    czz1       (m,)  <Sz_0 Sz_r> of the n_uc=1 solve at the largest maxm
    czz_a      (m,)  <Sz_0 Sz_r> of the n_uc=2 solve, from site 0 of the cell
    czz_b      (m,)  <Sz_1 Sz_{1+r}> of the n_uc=2 solve
    czz        (m,)  (czz_a + czz_b)/2, the correlator averaged over the cell
    e_uc2      ()    energy per site of the n_uc=2 solve
    params           JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import infinitechain

params = common.start(dict(
    model="infinite spin-1/2 Heisenberg chain, n_uc=1, J=1",
    gs_method="idmrg",
    itensor_version=int(os.environ.get("FIG_ITV", 3)),
    maxm=[int(x) for x in os.environ.get(
        "FIG_MAXM", "8 12 16 24 32 48 64 96").split()],
    maxiter=int(os.environ.get("FIG_MAXITER", 400)),
    etol=1e-10,
    rmax=int(os.environ.get("FIG_RMAX", 40)),
    correlator_cell="n_uc=2 at the largest maxm, averaged over the two "
                    "sites: the iDMRG state of the critical chain has a "
                    "period of two sites (a finite-chi dimerization)",
    seeds="none (iDMRG grows from a deterministic start)",
))
E_EXACT = 0.25 - np.log(2.0)


def chain(m):
    ic = infinitechain.Infinite_Spin_Chain(
        ["1/2"], itensor_version=params["itensor_version"])
    ic.set_hamiltonian(ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0]
                       + ic.SzC[0]*ic.SzR[0])
    ic.gs_method = params["gs_method"]
    ic.maxm = m
    ic.maxiter = params["maxiter"]
    ic.etol = params["etol"]
    return ic


energy, conv, times = [], [], []
for m in params["maxm"]:
    t0 = time.time()
    ic = chain(m)
    e = ic.gs_energy()
    energy.append(e)
    conv.append(bool(ic.converged))
    times.append(time.time() - t0)
    print("maxm=%3d  e=%.12f  e-e_exact=%.3e  converged=%s  (%.1f s)"
          % (m, e, e - E_EXACT, ic.converged, times[-1]), flush=True)

# correlators. iDMRG gives the critical chain a period of two sites at
# finite maxm (alternating bonds), so <Sz_0 Sz_1> of the n_uc=1 state is
# one of two bonds and 3<Sz_0 Sz_1> != e. On a two-site cell both bonds
# are measured, and their average gives back e.
r = np.arange(1, params["rmax"] + 1)
czz1 = np.array([ic.correlator("Sz", 0, "Sz", int(x)).real for x in r])
print("n_uc=1: 3<Sz0 Sz1> = %.10f  e = %.10f" % (3*czz1[0], energy[-1]))
t0 = time.time()
ic2 = infinitechain.Infinite_Spin_Chain(
    ["1/2", "1/2"], itensor_version=params["itensor_version"])
h = 0
for o in ("Sx", "Sy", "Sz"):
    oc, orr = getattr(ic2, o + "C"), getattr(ic2, o + "R")
    h = h + oc[0]*oc[1] + oc[1]*orr[0]
ic2.set_hamiltonian(h)
ic2.gs_method = params["gs_method"]
ic2.maxm = params["maxm"][-1]
ic2.maxiter = params["maxiter"]
ic2.etol = params["etol"]
e_uc2 = ic2.gs_energy()
czz_a = np.array([ic2.correlator("Sz", 0, "Sz", int(x)).real for x in r])
czz_b = np.array([ic2.correlator("Sz", 1, "Sz", int(x)).real for x in r])
czz = 0.5*(czz_a + czz_b)
print("n_uc=2: e=%.10f  bonds 3<SzSz> = %.10f %.10f  average %.10f  (%.1f s)"
      % (e_uc2, 3*czz_a[0], 3*czz_b[0], 3*czz[0], time.time() - t0))
for x, c1, ca, cb, c in zip(r, czz1, czz_a, czz_b, czz):
    print("r=%3d  n_uc=1 % .6e  site0 % .6e  site1 % .6e  avg % .6e"
          % (x, c1, ca, cb, c))

common.save(__file__, params, maxm=np.array(params["maxm"]),
            energy=np.array(energy), converged=np.array(conv),
            time_s=np.array(times), e_exact=np.array(E_EXACT),
            r=r, czz1=czz1, czz_a=czz_a, czz_b=czz_b, czz=czz,
            e_uc2=np.array(e_uc2))
