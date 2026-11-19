"""Should TASK expand its bonds by default on a chain with several baths?
One state of the spin dimer on two metals (dimer.py's model at Js = 0.3,
maxm 128), and on it: TASK with no expansion, with expand_bonds = 8 and
expand_rounds = 3, and with the two stronger settings 16/3 and 8/4, each
timed and with the weight it misses (taskdmrg.get_poles), and the
correction vector (CVM) at the same Lorentzian broadening as the
reference.

    python3 promote.py <outdir>

Writes <outdir>/promote.npz and prints one line per setting."""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr, taskdmrg
from dimer import LAMBDA, N, JK

ES = np.linspace(-0.6, 0.6, 121)
DELTA = 0.02
SETTINGS = [None, (8, 3), (16, 3), (8, 4)]


def run(outdir, order="sides", settings=SETTINGS, tag=""):
    fc = ysr.spin_chain_on_metal(2, N, site=[0, 1], Js=0.3, JK=JK,
                                 Lambda=LAMBDA, order=order)
    fc.maxm, fc.nsweeps, fc.cutoff = 128, 30, 1e-10
    e0 = float(np.real(fc.gs_energy()))
    p = fc.impurities[0]
    out = dict(es=ES, e0=e0)
    for s in settings:
        kw = dict(nkry=100)
        if s is not None:
            kw.update(expand_bonds=s[0], expand_rounds=s[1])
        t0 = time.time()
        d = ysr.get_stm_spectrum(fc, ES, delta=DELTA, imp=p, **kw)[1]
        t = time.time() - t0
        miss, dim = 0.0, 0
        for X in ysr.tunneling_operators(fc, imp=p):
            Xd = X.get_dagger()
            for name in ((Xd, X), (X, Xd)):
                r = taskdmrg.get_poles(fc, name=name, **kw)
                miss += abs(r["missing_weight"])
                dim = max(dim, r["dimension"])
        lab = "none" if s is None else "%d_%d" % s
        out["task_" + lab], out["time_" + lab] = d, t
        out["missing_" + lab], out["dim_" + lab] = miss, dim
        print(lab, "time %.0f s" % t, "missing %.3e" % miss, "dim", dim,
              "peak %.4f" % np.max(d), flush=True)
    t0 = time.time()
    es_c = ES[::6]
    c = ysr.get_stm_spectrum(fc, es_c, delta=DELTA, imp=p, submode="CVM")[1]
    out["es_cvm"], out["cvm"], out["time_cvm"] = es_c, c, time.time() - t0
    for s in settings:
        lab = "none" if s is None else "%d_%d" % s
        dd = np.max(np.abs(out["task_" + lab][::6] - c))/np.max(c)
        print(lab, "against CVM: %.3e of the peak" % dd, flush=True)
    np.savez(os.path.join(outdir, "promote%s.npz" % tag), **out)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[2] == "linear":
        # the tip spin listed first, so its chain comes right after the
        # spins: the vertex is then inside the tangent space
        run(sys.argv[1], order="linear", settings=[None, (8, 3)],
            tag="_linear")
    else:
        run(sys.argv[1])
