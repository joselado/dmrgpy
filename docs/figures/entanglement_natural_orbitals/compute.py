"""Figure 8, entanglement_natural_orbitals: natural-orbital occupations
n_alpha and correlation entropy S of the half-filled Hubbard chain as U
grows, with an ED check of S(U) on a 4-site chain.

Run with  bash docs/figures/run.sh entanglement_natural_orbitals

Arrays in data.npz:
    U        (k,)      Hubbard interactions, in units of t
    occ      (k,2n)    sorted occupations n_alpha of the n-site chain (DMRG)
    s        (k,)      correlation entropy of the n-site chain (DMRG)
    energy   (k,)      ground-state energy of the n-site chain (DMRG)
    s_small_ed   (k,)  correlation entropy of the 4-site chain, ED
    s_small_dmrg (k,)  correlation entropy of the 4-site chain, DMRG
    U_fine   (m,)      a finer U grid for the 4-site ED curve
    s_fine_ed (m,)     correlation entropy of the 4-site chain on U_fine, ED
    params             JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import fermionchain

n = int(os.environ.get("FIG_N", 12))
params = common.start(dict(
    model="Hubbard chain, open, H = -t sum_s (c+_is c_i+1s + h.c.) "
          "+ U sum Nup_i Ndn_i, t=1",
    n=n,
    n_small=4,
    U=[0, 1, 2, 4, 8, 16],
    U_fine=list(np.round(np.linspace(0, 16, 33), 6)),
    filling="half, set_conserved_sector(Nf=n)",
    maxm=100,
    nsweeps=24,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))


def chain(nsite, U):
    fc = fermionchain.Spinful_Fermionic_Chain(
        nsite, itensor_version=params["itensor_version"])
    h = 0
    for i in range(nsite-1):
        for c, cd in ((fc.Cup, fc.Cdagup), (fc.Cdn, fc.Cdagdn)):
            h = h - (cd[i]*c[i+1] + cd[i+1]*c[i])
    if U != 0:
        h = h + U*sum(fc.Nup[i]*fc.Ndn[i] for i in range(nsite))
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nsite)
    fc.maxm = params["maxm"]
    fc.nsweeps = params["nsweeps"]
    return fc


occ, s, energy, s_small_ed, s_small_dmrg = [], [], [], [], []
for U in params["U"]:
    t0 = time.time()
    fc = chain(n, U)
    e = fc.gs_energy()
    na = np.sort(fc.get_correlation_eigenvalues().real)[::-1]
    sc = fc.get_correlation_entropy()
    occ.append(na)
    s.append(sc)
    energy.append(e)
    f4 = chain(params["n_small"], U)
    sd = f4.get_correlation_entropy()
    se = f4.get_correlation_entropy(mode="ED")
    s_small_dmrg.append(sd)
    s_small_ed.append(se)
    print("U=%5.1f  E0=%.8f  S=%.6f  N=%.6f  | n=4: S_DMRG=%.8f S_ED=%.8f "
          "diff=%.1e  (%.1f s)" % (U, e, sc, np.sum(na), sd, se, sd-se,
                                   time.time()-t0), flush=True)

s_fine_ed = []
for U in params["U_fine"]:
    s_fine_ed.append(chain(params["n_small"], U).get_correlation_entropy(mode="ED"))

common.save(__file__, params, U=np.array(params["U"]), occ=np.array(occ),
            s=np.array(s), energy=np.array(energy),
            s_small_ed=np.array(s_small_ed),
            s_small_dmrg=np.array(s_small_dmrg),
            U_fine=np.array(params["U_fine"]), s_fine_ed=np.array(s_fine_ed))
