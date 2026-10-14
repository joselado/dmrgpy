"""Figure 12, mean-field_neel: the self-consistent mean-field decoupling
of the spin-1/2 Heisenberg chain as the fraction p of the exchange kept
exactly goes from 0 (pure mean field, a Neel state) to 1 (the exact
chain, no order).

Run with  bash docs/figures/run.sh mean-field_neel

References: at p=0 the classical Neel state, |<Sz_i>| = 1/2 and
E = -(n-1)/4; at p=1 the exact open chain, <Sz_i> = 0 and the energy of
a plain DMRG solve, close to the Bethe-ansatz estimate
E = n (1/4 - ln 2) + (pi - 1 - 2 ln 2)/4 (bulk plus surface energy).

Arrays in data.npz:
    p        (k,)     fractions of the exchange kept exactly
    mz       (k,n)    self-consistent <Sz_i>
    ms       (k,)     staggered magnetization of the central sites
    e_var    (k,)     <H> of the original Heisenberg chain in the
                      self-consistent state (a variational energy), total
    e_mf     (k,)     gs_energy() of the returned mean-field chain, total
    e_exact  ()       DMRG ground-state energy of the chain, total
    params            JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain, meanfield

n = int(os.environ.get("FIG_N", 24))
params = common.start(dict(
    model="spin-1/2 Heisenberg chain, open, J=1",
    n=n,
    p=[0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0],
    m0="staggered, [0, 0, 0.5 (-1)^i]",
    mix=0.5,
    maxerror=1e-6,
    maxite=500,
    central_sites=min(12, n//2),
    maxm=60,
    nsweeps=12,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))

sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=params["itensor_version"])
sc.maxm, sc.nsweeps = params["maxm"], params["nsweeps"]
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
e_exact = sc.gs_energy()
e_ba = n*(0.25-np.log(2)) + (np.pi-1-2*np.log(2))/4
print("exact chain: E0/n = %.6f, Bethe ansatz bulk+surface %.6f, "
      "bulk %.6f" % (e_exact/n, e_ba/n, 0.25-np.log(2)), flush=True)

m0 = [[0, 0, 0.5*(-1)**i] for i in range(n)]
c0 = (n - params["central_sites"])//2
centre = np.arange(c0, c0 + params["central_sites"])
mz, ms, e_var, e_mf = [], [], [], []
for p in params["p"]:
    t0 = time.time()
    mf = meanfield.spinchain_meanfield(sc, p=p, m0=m0, mix=params["mix"],
                                       maxerror=params["maxerror"],
                                       maxite=params["maxite"])
    m = np.array(mf.get_magnetization()).real      # 3 x n
    z = m[2]
    s = np.mean([(-1.0)**i*z[i] for i in centre])
    ev = mf.vev(h).real
    em = mf.gs_energy()
    mz.append(z)
    ms.append(s)
    e_var.append(ev)
    e_mf.append(em)
    print("p=%.2f  m_s=%.6f  <H>/n=%.6f  E_MF/n=%.6f  max|Sx|,|Sy|=%.1e  "
          "(%.1f s)" % (p, s, ev/n, em/n, np.max(np.abs(m[:2])),
                        time.time()-t0), flush=True)

common.save(__file__, params, p=np.array(params["p"]), mz=np.array(mz),
            ms=np.array(ms), e_var=np.array(e_var), e_mf=np.array(e_mf),
            e_exact=np.array(e_exact))
