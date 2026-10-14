"""Figure 1, ground-state_convergence: energy error and energy fluctuation
of the spin-1/2 Heisenberg chain against the bond dimension.

Run with  bash docs/figures/run.sh ground-state_convergence

Arrays in data.npz:
    maxm     (k,)  bond dimensions
    energy   (k,)  DMRG ground-state energy E0(maxm), total, in units of J
    de       (k,)  energy fluctuation dE = sqrt(<H^2>-<H>^2), total, in J
    e_ed     ()    ED ground-state energy, total, in J
    params         JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    n=16,
    maxm=[2, 4, 8, 12, 16, 24, 32, 48, 64, 96],
    nsweeps=20,
    cutoff=1e-14,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]


def chain():
    sc = spinchain.Spin_Chain(["S=1/2"] * n,
                              itensor_version=params["itensor_version"])
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    return sc


e_ed = chain().gs_energy(mode="ED")
print("ED  E0 = %.12f" % e_ed, flush=True)

energy, de = [], []
for m in params["maxm"]:
    t0 = time.time()
    sc = chain()          # a fresh chain per maxm, so nothing is reused
    sc.maxm = m
    sc.nsweeps = params["nsweeps"]
    sc.cutoff = params["cutoff"]
    e = sc.gs_energy()
    d = sc.gs_energy_fluctuation()
    energy.append(e)
    de.append(d)
    print("maxm=%3d  E0=%.12f  E0-E_ED=%.3e  dE=%.3e  (%.1f s)"
          % (m, e, e - e_ed, d, time.time() - t0), flush=True)

common.save(__file__, params, maxm=np.array(params["maxm"]),
            energy=np.array(energy), de=np.array(de), e_ed=np.array(e_ed))
