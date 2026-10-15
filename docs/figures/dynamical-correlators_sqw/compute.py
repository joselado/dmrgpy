# Dynamical spin structure factor S^zz(q,w) of the open spin-1/2 Heisenberg
# chain, one KPM run per momentum. On an open chain the momentum states are
# standing waves, S^z_q = sqrt(2/(n+1)) sum_i sin(q(i+1)) S^z_i, q = pi k/(n+1).
import numpy as np
from dmrgpy import spinchain, kpmdmrg

n = 32
delta = 0.15
es = np.linspace(0.0, 4.0, 400)
qs = np.pi*np.arange(1, n+1)/(n+1)

sc = spinchain.Spin_Chain(["S=1/2"]*n)
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
sc.maxm, sc.nsweeps, sc.cutoff = 80, 14, 1e-12
sc.kpmmaxm = 40  # kpmmaxm=80 changes S(q,w) at q close to pi by less than 1%
sc.gs_energy()


def spectrum(x, mus, emin, emax, scale, npol, d):
    _, y = kpmdmrg.dynamical_correlator_from_moments(
        mus, emin, emax, scale, npol, x, kernel="jackson", delta=d,
        origin=sc.e0)
    return y.real


sqw, sumrule, sq = [], [], []
for q in qs:
    Sq = sum(np.sqrt(2/(n+1))*np.sin(q*(i+1))*sc.Sz[i] for i in range(n))
    m = sc.get_dynamical_correlator_moments(name=(Sq, Sq), delta=delta)
    sqw.append(spectrum(es, *m))
    # sum rule: the integral over the whole spectrum is <S^z_q S^z_q>
    ew = np.linspace(-1.0, m[2] - sc.e0 + 1.0, 20000)
    sumrule.append(np.trapezoid(spectrum(ew, *m), ew))
    sq.append(sc.vev(Sq*Sq).real)

np.savez("data.npz", qs=qs, es=es, sqw=sqw, sumrule=sumrule, sq=sq)
