"""A chain of ns spins 1/2 with every spin on its own metal, the dI/dV under
one spin by plain TASK in the linear order with that spin listed first
(its chain right after the spins, the only layout whose tangent space
holds the tip vertex without bond expansion). One chain is solved per tip.

    python3 chain_linear.py <task> <outdir>

The tasks are the rows of TASKS (ns, tip, bond dimension). Every task
writes <outdir>/chain_linear_<task>.npz with the bias, dI/dV, the
ground-state energy, the weight of the vertices TASK misses, the setting
and the wall time."""
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr
from chain import LAMBDA, N, JK, ES, DELTA

TASKS = [dict(ns=ns, tip=tip, maxm=maxm)
         for ns in (3, 4) for tip in (0, 1) for maxm in (128, 256)]


def run(k, outdir):
    s = TASKS[k]
    t0 = time.time()
    site = [s["tip"]] + [i for i in range(s["ns"]) if i != s["tip"]]
    fc = ysr.spin_chain_on_metal(s["ns"], N, site=site, Js=0.3, JK=JK,
                                 Lambda=LAMBDA, order="linear")
    fc.maxm, fc.nsweeps, fc.cutoff = s["maxm"], 30, 1e-10
    e0 = float(np.real(fc.gs_energy()))
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        d = ysr.get_stm_spectrum(fc, ES, delta=DELTA, imp=fc.impurities[0],
                                 nkry=200)[1]
    miss = [str(w.message)[:160] for w in rec if "miss" in str(w.message)]
    np.savez(os.path.join(outdir, "chain_linear_%d.npz" % k), es=ES,
             didv=d, e0=e0, missing=json.dumps(miss), setting=json.dumps(s),
             time=time.time() - t0)
    print(k, s, "E0 = %.10f" % e0, "t = %.0f s" % (time.time() - t0))
    for m in miss:
        print("  ", m)


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
