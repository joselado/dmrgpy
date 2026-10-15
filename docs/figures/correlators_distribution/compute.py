"""Figure 9, correlators_distribution: the distribution P(m) of the
magnetization per site m = (1/n) sum_i Sz_i in the ground state of the
transverse-field Ising chain, in the ordered phase, at the critical field
and in the paramagnet.

Run with  bash docs/figures/run.sh correlators_distribution

The check is done on the small chain: DMRG and ED with the same Chebyshev
expansion (same window, same number of moments, same Jackson kernel), so
the two curves must agree point by point. The large chain is DMRG only.

Arrays in data.npz:
    h          (k,)    transverse fields
    m_small    (p,)    grid of m for the small chain
    p_small    (k,p)   P(m), small chain, DMRG
    p_small_ed (k,p)   P(m), small chain, ED
    m_large    (q,)    grid of m for the large chain
    p_large    (k,q)   P(m), large chain, DMRG
    mz_large   (k,)    <m> of the large chain from vev
    norm_large (k,)    integral of P(m) dm, large chain
    params             JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import spinchain

params = common.start(dict(
    model="H = -sum SzSz - h sum Sx - hz sum Sz, open, critical at h=1/2",
    h=[0.2, 0.5, 1.0],
    hz=1e-4,
    n_small=int(os.environ.get("FIG_NS", 12)),
    n_large=int(os.environ.get("FIG_NL", 24)),
    maxm=60,
    nsweeps=12,
    kpmmaxm=60,
    scale=0.625,
    npol=96,
    kpm_window="X/scale, so the spectrum [-1/2,1/2] of m maps to "
               "[-0.8,0.8]; npol Chebyshev moments, Jackson kernel, on both "
               "solvers (DMRG: delta=3*scale/npol; ED: delta=4*scale/npol)",
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
scale, npol = params["scale"], params["npol"]
# int() inside both solvers rounds down, so ask for half a moment more;
# the DMRG route returns int(3*scale/delta)+2 moments, ED 4*int(scale/delta)
delta = 3*scale/(npol-2+0.5)
delta_ed = scale/(npol/4+0.5)


def chain(n, h):
    sc = spinchain.Spin_Chain(["S=1/2"]*n,
                              itensor_version=params["itensor_version"])
    H = -sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    H = H - h*sum(sc.Sx[i] for i in range(n))
    H = H - params["hz"]*sum(sc.Sz[i] for i in range(n))
    sc.set_hamiltonian(H)
    sc.maxm = params["maxm"]
    sc.kpmmaxm = params["kpmmaxm"]
    sc.nsweeps = params["nsweeps"]
    return sc


p_small, p_small_ed, p_large, mz_large, norm_large = [], [], [], [], []
for h in params["h"]:
    t0 = time.time()
    n = params["n_small"]
    sc = chain(n, h)
    M = sum(sc.Sz[i] for i in range(n))*(1.0/n)
    grid, p = sc.get_distribution(X=M, delta=delta, scale=scale)
    xe, pe = sc.get_distribution(X=M, delta=delta_ed, scale=scale,
                                 xs=grid, mode="ED")
    p, pe = p.real, pe.real
    p_small.append(p)
    p_small_ed.append(pe)
    print("h=%.2f  n=%d  max|P_DMRG-P_ED|=%.2e  max P=%.3f  int P_DMRG=%.6f"
          "  int P_ED=%.6f  (%.1f s)"
          % (h, n, np.max(np.abs(p-pe)), np.max(pe), np.trapezoid(p, grid),
             np.trapezoid(pe, grid), time.time()-t0), flush=True)
    t0 = time.time()
    n = params["n_large"]
    sc = chain(n, h)
    M = sum(sc.Sz[i] for i in range(n))*(1.0/n)
    grid_l, p = sc.get_distribution(X=M, delta=delta, scale=scale)
    p = p.real
    mz = sc.vev(M).real
    nrm = np.trapezoid(p, grid_l)
    p_large.append(p)
    mz_large.append(mz)
    norm_large.append(nrm)
    print("h=%.2f  n=%d  <m>=%.5f  int P=%.6f  int m P=%.5f  peak at m=%.4f"
          "  (%.1f s)" % (h, n, mz, nrm, np.trapezoid(grid_l*p, grid_l),
                          grid_l[np.argmax(p)], time.time()-t0), flush=True)

common.save(__file__, params, h=np.array(params["h"]), m_small=grid,
            p_small=np.array(p_small), p_small_ed=np.array(p_small_ed),
            m_large=grid_l, p_large=np.array(p_large),
            mz_large=np.array(mz_large), norm_large=np.array(norm_large))
