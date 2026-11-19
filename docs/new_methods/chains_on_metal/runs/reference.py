"""The reference dI/dV of the spin dimer on two metals (dimer.py's model at
Js = 0.3, the tip on the first spin) from the variational correction
vector (cvm_solver="variational", dynamical DMRG), at the biases and the
broadening of promote.py. The default conjugate-gradient CVM did not
converge on these vertices (11% of the peak off ED on 8 sites, where the
variational one is 6e-9 off).

    python3 reference.py <task> <outdir>

A task is one order, one bond dimension and one chunk of three biases;
every task writes <outdir>/reference_<task>.npz."""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr
from dimer import LAMBDA, N, JK

ES = np.linspace(-0.6, 0.6, 121)[::6] # the 21 biases of promote.py's CVM
DELTA = 0.02
CHUNKS = [ES[i:i + 3] for i in range(0, len(ES), 3)]
TASKS = [dict(order=o, maxm=m, chunk=c) for (o, m) in
         (("sides", 128), ("linear", 128), ("sides", 256))
         for c in range(len(CHUNKS))]


def run(k, outdir):
    s = TASKS[k]
    t0 = time.time()
    fc = ysr.spin_chain_on_metal(2, N, site=[0, 1], Js=0.3, JK=JK,
                                 Lambda=LAMBDA, order=s["order"])
    fc.maxm, fc.nsweeps, fc.cutoff = s["maxm"], 30, 1e-10
    fc.cvm_solver, fc.cvm_maxm = "variational", s["maxm"]
    e0 = float(np.real(fc.gs_energy()))
    es = CHUNKS[s["chunk"]]
    d = ysr.get_stm_spectrum(fc, es, delta=DELTA, imp=fc.impurities[0],
                             submode="CVM")[1]
    np.savez(os.path.join(outdir, "reference_%d.npz" % k), es=es, didv=d,
             e0=e0, setting=json.dumps(s), time=time.time() - t0)
    print(k, s, "E0 = %.10f" % e0, "t = %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
