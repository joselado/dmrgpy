"""The spin dimer on two metals of dimer.py, with plain TASK in the linear
order and the tip spin listed first: the Js sweep and the bond dimension.

    python3 dimer_linear.py <task> <outdir>

Every task writes <outdir>/dimer_linear_<task>.npz with the bias, the
dI/dV at the broadenings of dimer.py, the ground-state energy, the TASK
missing-weight warnings, the setting and the wall time."""
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr
from dimer import LAMBDA, N, JK, ES, DELTAS

TASKS = ([dict(Js=0.3, maxm=m, nkry=200) for m in (64, 128, 256)]
         + [dict(Js=js, maxm=128, nkry=200)
            for js in (0.0, 0.003, 0.01, 0.03, 0.1, 0.2, 0.5)])


def run(k, outdir):
    s = dict(TASKS[k])
    t0 = time.time()
    fc = ysr.spin_chain_on_metal(2, N, site=[0, 1], Js=s.pop("Js"), JK=JK,
                                 Lambda=LAMBDA, order="linear")
    fc.maxm, fc.nsweeps, fc.cutoff = s.pop("maxm"), 30, 1e-10
    e0 = float(np.real(fc.gs_energy()))
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        d = [ysr.get_stm_spectrum(fc, ES, delta=dl, imp=fc.impurities[0],
                                  **s)[1] for dl in DELTAS]
    miss = [str(w.message)[:160] for w in rec if "miss" in str(w.message)]
    np.savez(os.path.join(outdir, "dimer_linear_%d.npz" % k), es=ES,
             didv=d, deltas=DELTAS, e0=e0, missing=json.dumps(miss),
             setting=json.dumps(TASKS[k]), time=time.time() - t0)
    print(k, TASKS[k], "E0 = %.10f" % e0, "t = %.0f s" % (time.time() - t0))
    for m in miss:
        print("  ", m)


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
