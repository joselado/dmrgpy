"""A spin dimer on two metals, each spin with its own Wilson chain: the
dI/dV under the first spin, for one setting of the solver per call.

    python3 dimer.py <task> <outdir>

The tasks are the rows of TASKS: the first ones measure the convergence
of TASK (bond dimension, Krylov size, bond expansion) against the
correction vector (CVM, the same Lorentzian broadening) at
Js = 0.3, the others sweep Js at the converged setting. Every task
writes <outdir>/dimer_<task>.npz with the bias, dI/dV, the ground-state
energy, the setting and the wall time."""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr

LAMBDA, N = 2.0, 21 # Wilson chains reach about 2^-10 of the band
RHOJ = 0.15 # rho J in the continuum, JK = 2 D rho J on the chain times A
JK = 2*RHOJ*ysr.lambda_factor(LAMBDA)
ES = np.linspace(-0.6, 0.6, 241)
DELTAS = (0.01, 0.02) # Lorentzian half-widths of the TASK curves
ES_CVM, DELTA_CVM = ES[::10], 0.02 # the correction-vector check
X = dict(expand_bonds=8, expand_rounds=3)
TASKS = [ # convergence at Js = 0.3
    dict(Js=0.3, maxm=64, nkry=100, **X),
    dict(Js=0.3, maxm=128, nkry=100, **X),
    dict(Js=0.3, maxm=128, nkry=200, **X),
    dict(Js=0.3, maxm=256, nkry=200, **X),
    dict(Js=0.3, maxm=128, nkry=100), # no bond expansion
    dict(Js=0.3, maxm=128, submode="CVM"), # Lorentzian, as TASK
    dict(Js=0.3, maxm=256, submode="CVM"),
] + [dict(Js=js, maxm=128, nkry=200, **X) # the Js sweep
     for js in (0.0, 0.003, 0.01, 0.03, 0.1, 0.2, 0.5)]


def run(k, outdir):
    s = dict(TASKS[k])
    t0 = time.time()
    fc = ysr.spin_chain_on_metal(2, N, site=[0, 1], Js=s.pop("Js"), JK=JK,
                                 Lambda=LAMBDA, order="sides")
    fc.maxm, fc.nsweeps, fc.cutoff = s.pop("maxm"), 30, 1e-10
    e0 = float(np.real(fc.gs_energy()))
    if s.get("submode") == "CVM":
        es, deltas = ES_CVM, (DELTA_CVM,)
    else:
        es, deltas = ES, DELTAS
    d = [ysr.get_stm_spectrum(fc, es, delta=dl, imp=fc.impurities[0], **s)[1]
         for dl in deltas] # TASK builds its spectrum once and caches it
    np.savez(os.path.join(outdir, "dimer_%d.npz" % k), es=es, didv=d,
             deltas=deltas, e0=e0, setting=json.dumps(TASKS[k]),
             time=time.time() - t0)
    print(k, TASKS[k], "E0 = %.10f" % e0, "t = %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
