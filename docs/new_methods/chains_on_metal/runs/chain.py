"""A chain of ns spins 1/2 with every spin on its own metal: the dI/dV
under the end spin and under the next one, at two bond dimensions.

    python3 chain.py <task> <outdir>

The tasks are the rows of TASKS (ns, order, bond dimension). Every task
writes <outdir>/chain_<task>.npz with the bias, the two dI/dV curves,
the ground-state energy, the largest bond entropy, the setting and the
wall time."""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr

LAMBDA, N = 2.0, 15
RHOJ = 0.15
JK = 2*RHOJ*ysr.lambda_factor(LAMBDA)
ES = np.linspace(-1.0, 1.0, 201)
DELTA = 0.02
X = dict(expand_bonds=8, expand_rounds=3, nkry=200)
TASKS = [dict(ns=ns, order=order, maxm=maxm)
         for ns in (3, 4) for order in ("sides", "comb")
         for maxm in (128, 256)]


def run(k, outdir):
    s = TASKS[k]
    t0 = time.time()
    fc = ysr.spin_chain_on_metal(s["ns"], N, site=list(range(s["ns"])),
                                 Js=0.3, JK=JK, Lambda=LAMBDA,
                                 order=s["order"])
    fc.maxm, fc.nsweeps, fc.cutoff = s["maxm"], 30, 1e-10
    e0 = float(np.real(fc.gs_energy()))
    wf = fc.get_gs()
    smax = max(float(wf.get_bond_entropy(b)) for b in range(fc.ns - 1))
    d = [ysr.get_stm_spectrum(fc, ES, delta=DELTA, imp=fc.spin_sites[i],
                              **X)[1] for i in (0, 1)]
    np.savez(os.path.join(outdir, "chain_%d.npz" % k), es=ES, didv=d, e0=e0,
             smax=smax, setting=json.dumps(s), time=time.time() - t0)
    print(k, s, "E0 = %.10f" % e0, "Smax = %.3f" % smax,
          "t = %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
