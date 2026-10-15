"""Figure 22, post-processing_extrapolation: energy per site of open
spin-1/2 Heisenberg chains against 1/L, extrapolated to L -> infinity
with extrapolate.size_extrapolation.

Run with  bash docs/figures/run.sh post-processing_extrapolation

Arrays in data.npz:
    L         (k,)  chain lengths
    energy    (k,)  DMRG ground-state energy, total, in units of J
    time_s    (k,)  wall time of each solve, in s
    e_inf_n1  ()    size_extrapolation(L, E/L, n=1)
    e_inf_n2  ()    size_extrapolation(L, E/L, n=2)
    e_exact   ()    Bethe-ansatz energy per site 1/4 - ln 2
    params          JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, extrapolate

params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    L=[int(x) for x in os.environ.get(
        "FIG_L", "16 24 32 48 64 96 128").split()],
    maxm=int(os.environ.get("FIG_MAXM", 100)),
    nsweeps=16,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
E_EXACT = 0.25 - np.log(2.0)

energy, times = [], []
for L in params["L"]:
    t0 = time.time()
    sc = spinchain.Spin_Chain(["S=1/2"] * L,
                              itensor_version=params["itensor_version"])
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(L-1))
    sc.set_hamiltonian(h)
    sc.maxm = params["maxm"]
    sc.nsweeps = params["nsweeps"]
    e = sc.gs_energy()
    energy.append(e)
    times.append(time.time() - t0)
    print("L=%4d  E=%.10f  E/L=%.10f  (%.1f s)" % (L, e, e/L, times[-1]),
          flush=True)

L = np.array(params["L"])
energy = np.array(energy)
e1 = extrapolate.size_extrapolation(L, energy/L, n=1)
e2 = extrapolate.size_extrapolation(L, energy/L, n=2)
print("n=1: e_inf=%.8f  error %.2e" % (e1, e1 - E_EXACT))
print("n=2: e_inf=%.8f  error %.2e" % (e2, e2 - E_EXACT))
print("exact      %.8f" % E_EXACT)

common.save(__file__, params, L=L, energy=energy, time_s=np.array(times),
            e_inf_n1=np.array(e1), e_inf_n2=np.array(e2),
            e_exact=np.array(E_EXACT))
