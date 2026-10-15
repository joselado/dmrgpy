"""Figure 11, fidelity_ising: fidelity susceptibility per site of the
open transverse-field Ising chain H = sum sz_i sz_i+1 + lambda sum sx_i
(Pauli matrices) across its critical point lambda = 1.

Run with  bash docs/figures/run.sh fidelity_ising

References: mode="ED" at the smallest length, and the bulk result of the
infinite chain, chi/n = 1/(16(1-l^2)) for l<1 and 1/(16 l^2 (l^2-1))
for l>1, which the curves approach away from l=1.

Arrays in data.npz:
    lam      (k,)    lambda grid
    ns       (p,)    chain lengths
    chi      (p,k)   fidelity susceptibility chi (total, not per site), DMRG
    lam_ed   (k,)    lambda grid of the ED reference (same as lam)
    chi_ed   (k,)    chi of the smallest chain, ED
    params           JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, fidelity

ns_env = os.environ.get("FIG_NS")
params = common.start(dict(
    model="H = h0 + lambda h1, h0 = 4 sum SzSz, h1 = 2 sum Sx, open",
    ns=[int(a) for a in ns_env.split(",")] if ns_env else [12, 16, 24, 32],
    lam=list(np.round(np.linspace(0.8, 1.4, int(os.environ.get("FIG_NL", 49))), 8)),
    dl=0.02,
    fmode="derivative",
    tol=float(os.environ.get("FIG_TOL", 1e-4)),
    maxm=60,
    nsweeps=20,
    note="lambda starts at 0.8: below about 0.75 the Z2 doublet of the "
         "ordered open chain splits by less than tol=1e-4 at n=32, and the "
         "derivative estimator then returns chi summed over both states of "
         "the doublet; with a smaller tol DMRG lands on an arbitrary "
         "superposition of the doublet and the overlaps are meaningless",
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))


def chain(n):
    sc = spinchain.Spin_Chain(["S=1/2"]*n,
                              itensor_version=params["itensor_version"])
    sc.maxm, sc.nsweeps = params["maxm"], params["nsweeps"]
    h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    h1 = 2*sum(sc.Sx[i] for i in range(n))
    sc.set_hamiltonian(h0 + h1)
    return sc, h0, h1


def bulk(l):
    l = np.asarray(l)
    return np.where(l < 1, 1/(16*(1-l**2)), 1/(16*l**2*(l**2-1)))


chi = []
for n in params["ns"]:
    t0 = time.time()
    sc, h0, h1 = chain(n)
    row = [fidelity.get_fidelity(sc, h0, h1, l, dl=params["dl"],
                                 delta=params["tol"], fmode=params["fmode"])
           for l in params["lam"]]
    chi.append(row)
    k = int(np.argmax(row))
    print("n=%d  peak chi/n=%.4f at lambda=%.3f;  chi/n at %.2f: %.5f, "
          "at 1.4: %.5f (bulk %.5f)  (%.1f s)"
          % (n, row[k]/n, params["lam"][k], params["lam"][0], row[0]/n,
             row[-1]/n, bulk(1.4), time.time()-t0), flush=True)

# the edge contribution cancels in a difference of two lengths
c = np.array(chi)
dens = (c[-1] - c[-2])/(params["ns"][-1] - params["ns"][-2])
for j in (0, len(params["lam"])//2, -1):
    l = params["lam"][j]
    print("bulk density from n=%d,%d at lambda=%.3f: %.5f, exact infinite "
          "chain %.5f" % (params["ns"][-2], params["ns"][-1], l, dens[j],
                          bulk(l)), flush=True)

n = params["ns"][0]
sc, h0, h1 = chain(n)
chi_ed = np.array([fidelity.get_fidelity(sc, h0, h1, l, dl=params["dl"],
                                         delta=params["tol"],
                                         fmode=params["fmode"], mode="ED")
                   for l in params["lam"]])
print("ED n=%d: max|chi_DMRG-chi_ED|/max chi_ED = %.2e" %
      (n, np.max(np.abs(np.array(chi[0])-chi_ed))/np.max(chi_ed)), flush=True)

common.save(__file__, params, lam=np.array(params["lam"]),
            ns=np.array(params["ns"]), chi=np.array(chi),
            lam_ed=np.array(params["lam"]), chi_ed=chi_ed)
