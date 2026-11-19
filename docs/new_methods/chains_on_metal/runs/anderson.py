"""An Anderson level at the end of a chain of spins 1/2, on a metal: dI/dV
on the level by the NRG route (DMRG head, NRG tail, TASK above the
hand-over), averaged over twists of the mesh, for one setting per call.
The environment sets nz (NZ, 4), the bond dimension of the head (MAXM,
200), the depth of the hand-over (M, 4) and the states kept (KEEP, 800);
seven spins need a hand-over basis of 300 to 400 states (MAXM=512, M=5,
docs/new_methods/kondo_nrg_mps/README.md).

    python3 anderson.py <task> <outdir>

The tasks are the rows of TASKS: one level alone at particle-hole
symmetry (the Friedel sum rule, pi Gamma rho_d(0) = 1 per spin), the
level at the end of 3 and 7 spins for three values of ed, and the Kondo
spin of spin_chain_on_metal at the Schrieffer-Wolff coupling of the
symmetric level, for comparison. Every task writes
<outdir>/anderson_<task>.npz with the bias, dI/dV, the setting and the
wall time."""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("DMRGPY_SRC", "../../../../src"))
import numpy as np
from dmrgpy import ysr

LAMBDA, N_HEAD, D = 4.0, 12, 1.0
U, GAMMA, JS = 0.5, 0.05, 0.01 # Js below T_K and above it are both of use
A = ysr.lambda_factor(LAMBDA)
TD = np.sqrt(2*GAMMA/np.pi) # Gamma = pi rho tD^2, rho = 1/(2D)
JK = 8*TD**2/U # Schrieffer-Wolff at ed = -U/2
ES = np.concatenate([-np.logspace(0, -6, 121), [0.0], np.logspace(-6, 0, 121)])
TASKS = [dict(ns=1, ed=-U/2, kind="anderson")]
TASKS += [dict(ns=ns, ed=ed, kind="anderson")
          for ns in (3, 7) for ed in (-U/2, -U/3, -U/6)]
TASKS += [dict(ns=ns, kind="kondo") for ns in (1, 3, 7)]


def run(k, outdir, nz=int(os.environ.get("NZ", 4)),
        maxm=int(os.environ.get("MAXM", 200)), m=int(os.environ.get("M", 4)),
        keep=int(os.environ.get("KEEP", 800))):
    s = TASKS[k]
    t0 = time.time()
    if s["kind"] == "anderson":
        fc = ysr.anderson_spin_chain(s["ns"], N_HEAD, Js=JS, ed=s["ed"], U=U,
                                     tD=TD*np.sqrt(A), Lambda=LAMBDA, D=D,
                                     itensor_version=3)
    else:
        fc = ysr.spin_chain_on_metal(s["ns"], N_HEAD, Js=JS, JK=JK*A,
                                     Lambda=LAMBDA, D=D, itensor_version=3)
    fc.maxm, fc.nsweeps, fc.cutoff = maxm, 30, 1e-12
    _, didv = ysr.get_stm_spectrum(fc, ES, submode="NRG", nz=nz, m=m,
                                   keep=keep)
    np.savez(os.path.join(outdir, "anderson_%d.npz" % k), es=ES, didv=didv,
             setting=json.dumps(dict(s, nz=nz, maxm=maxm, m=m, keep=keep)),
             gamma=GAMMA, U=U, time=time.time() - t0)
    print(k, s, "dI/dV(0) = %.4f" % didv[len(ES)//2],
          "t = %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    run(int(sys.argv[1]), sys.argv[2])
