# Distribution P(m) of the magnetization per site m in the ground state of the
# transverse-field Ising chain H = -sum SzSz - h sum Sx - hz sum Sz, critical
# at h=1/2, from a Chebyshev (KPM) expansion; ED checks DMRG on the small chain.
import numpy as np
from dmrgpy import spinchain

hs = [0.2, 0.5, 1.0]
hz = 1e-4  # tiny longitudinal field that picks one of the two ordered states
n_small, n_large = 12, 24
scale, npol = 0.625, 96  # m in [-1/2,1/2] maps to [-0.8,0.8], 96 moments
# The same npol moments on both solvers: DMRG keeps int(3*scale/delta)+2 and
# ED 4*int(scale/delta), and half a moment more survives the int() rounding.
delta = 3*scale/(npol - 2 + 0.5)
delta_ed = scale/(npol/4 + 0.5)


def chain(n, h):
    """The chain and its magnetization per site m."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    sc.set_hamiltonian(-sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
                       - h*sum(sc.Sx[i] for i in range(n))
                       - hz*sum(sc.Sz[i] for i in range(n)))
    sc.maxm, sc.kpmmaxm, sc.nsweeps = 60, 60, 12
    return sc, sum(sc.Sz[i] for i in range(n))*(1.0/n)


p_small, p_small_ed, p_large = [], [], []
for h in hs:
    sc, m = chain(n_small, h)
    m_small, p = sc.get_distribution(X=m, delta=delta, scale=scale)
    p_ed = sc.get_distribution(X=m, delta=delta_ed, scale=scale, xs=m_small,
                               mode="ED")[1]
    p_small.append(p.real)
    p_small_ed.append(p_ed.real)
    sc, m = chain(n_large, h)
    m_large, p = sc.get_distribution(X=m, delta=delta, scale=scale)
    p_large.append(p.real)

np.savez("data.npz", h=hs, n_small=n_small, n_large=n_large, m_small=m_small,
         p_small=p_small, p_small_ed=p_small_ed, m_large=m_large,
         p_large=p_large)
