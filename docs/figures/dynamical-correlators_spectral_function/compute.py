"""Figure 14, dynamical-correlators_spectral_function: the Mott gap of the
half-filled Hubbard chain opening with U.

Left panel: the local spectral function A_ii(omega) of a bulk site, spin
up, from the two KPM correlators (C_up, Cdag_up) (addition) and
(Cdag_up, C_up) (removal), on the particle-hole symmetric Hamiltonian

    H = -t sum_{i,s} (c^dag_{is} c_{i+1,s} + h.c.)
        + U sum_i (n_{i,up} - 1/2)(n_{i,dn} - 1/2),

whose global ground state is the half-filled singlet with mu = 0, so
omega is already omega - mu. At U=0 the exact free-fermion spectrum is
built from the exact Chebyshev moments on the same rescaling and the same
kernel, so the comparison checks the DMRG moments themselves.

Right panel: the charge gap E0(N+1)+E0(N-1)-2E0(N) and the chemical
potential mu = (E0(N+1)-E0(N-1))/2 from three conserved-sector ground
states of H' = -t sum (...) + U sum n_up n_dn (no 1/2 shifts), against the
Lieb-Wu gap of the infinite chain and mu = U/2.

Run with  bash docs/figures/run.sh dynamical-correlators_spectral_function

Arrays in data.npz:
    Us_spec    (nu,)       U of the spectra, in t
    es         (ne,)       omega - mu, in t
    A          (nu, ne)    A_ii(omega - mu), Jackson kernel, real part
    sumrule    (nu,)       int A domega over the full spectrum (exact: 1)
    A_free     (ne,)       exact free-fermion A_ii, same moments count, kernel
    A_free_lor (ne,)       exact free-fermion A_ii with Lorentzian lines
    mom_dev    ()          max |mu_m(DMRG) - mu_m(exact)| at U=0, both sides
    Us_gap     (ng,)       U of the gap scan, in t
    ns_gap     (nn,)       chain lengths of the gap scan
    gap        (nn, ng)    charge gap of the finite chains, in t
    mu         (nn, ng)    chemical potential of the finite chains, in t
    e_sectors  (nn, ng, 3) E0 of the sectors N-1, N, N+1, in t
    Ul         (nl,)       fine U grid of the Lieb-Wu curve
    gap_lw     (nl,)       Lieb-Wu charge gap of the infinite chain, in t
    gap_free   ()          free-fermion gap of the n-site chain, in t
    params                 JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import fermionchain, kpmdmrg

params = common.start(dict(
    model="Hubbard chain, open, hopping -t, t=1, half filling",
    n=int(os.environ.get("FIG_N", 12)),
    site=None,
    spin="up",
    Us_spec=[0.0, 2.0, 4.0, 8.0],
    Us_gap=[0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0],
    ns_gap=None,
    maxm=80,
    nsweeps=12,
    cutoff=1e-12,
    kpmmaxm=80,
    delta=0.25,
    kernel="jackson",
    es=[-8.0, 8.0, 800],
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
n = params["n"]
site = n // 2
params["site"] = site
params["ns_gap"] = [n, 2*n]
delta = params["delta"]
es = np.linspace(*params["es"][:2], params["es"][2])


def hubbard(u, shifted, n=n):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(
        n, itensor_version=params["itensor_version"])
    h = 0
    for i in range(n-1):
        for cd, c in ((fc.Cdagup, fc.Cup), (fc.Cdagdn, fc.Cdn)):
            h = h - cd[i]*c[i+1] - cd[i+1]*c[i]
    for i in range(n):
        if shifted:
            h = h + u*(fc.Nup[i] - 0.5)*(fc.Ndn[i] - 0.5)
        else:
            h = h + u*fc.Nup[i]*fc.Ndn[i]
    fc.set_hamiltonian(h)
    fc.maxm = params["maxm"]
    fc.nsweeps = params["nsweeps"]
    fc.cutoff = params["cutoff"]
    fc.kpmmaxm = params["kpmmaxm"]
    return fc


# free fermions on the same open chain: levels and weights at the site
hop = -(np.eye(n, k=1) + np.eye(n, k=-1))
eps, phi = np.linalg.eigh(hop)
w = phi[site, :]**2
gap_free = eps[n//2] - eps[n//2-1]
part = eps > 0     # empty levels at half filling (no zero level for even n)


def moments(fc, A, B):
    """KPM moments of (A,B) and the spectrum on es and on a wide grid."""
    mus, emin, emax, scale, npol, d = fc.get_dynamical_correlator_moments(
        name=(A, B), delta=delta)
    return dict(mus=np.asarray(mus), emin=emin, emax=emax, scale=scale,
                npol=npol, d=d, e0=fc.e0)


def profile(m, x, mus=None):
    _, y = kpmdmrg.dynamical_correlator_from_moments(
        m["mus"] if mus is None else mus, m["emin"], m["emax"], m["scale"],
        m["npol"], x, kernel=params["kernel"], delta=m["d"], origin=m["e0"])
    return y.real


def exact_moments(m, poles, weights):
    """Chebyshev moments sum_k w_k T_j(x_k) on the rescaling of m."""
    x = (poles + m["e0"] - 0.5*(m["emin"] + m["emax"]))*m["scale"]
    j = np.arange(len(m["mus"]))
    return np.sum(weights[None, :]*np.cos(j[:, None]*np.arccos(x[None, :])),
                  axis=1)


A_spec, sumrule = [], []
A_free = A_free_lor = None
mom_dev = np.nan
for u in params["Us_spec"]:
    t0 = time.time()
    fc = hubbard(u, shifted=True)
    e0 = fc.gs_energy()
    Nt = sum(fc.vev(fc.Nup[i] + fc.Ndn[i]).real for i in range(n))
    mp = moments(fc, fc.Cup[site], fc.Cdagup[site])     # addition, omega>0
    mh = moments(fc, fc.Cdagup[site], fc.Cup[site])     # removal, omega<0
    a = profile(mp, es) + profile(mh, -es)
    wide_p = np.linspace(-1.0, mp["emax"] - mp["e0"] + 1.0, 20000)
    wide_h = np.linspace(-1.0, mh["emax"] - mh["e0"] + 1.0, 20000)
    s = np.trapezoid(profile(mp, wide_p), wide_p) + \
        np.trapezoid(profile(mh, wide_h), wide_h)
    A_spec.append(a)
    sumrule.append(s)
    print("U=%.1f  E0=%.8f  <N>=%.6f  sum rule %.6f  npol %d/%d  (%.1f s)"
          % (u, e0, Nt, s, mp["npol"], mh["npol"], time.time()-t0), flush=True)
    if u == 0.0:
        exp_p = exact_moments(mp, eps[part], w[part])
        exp_h = exact_moments(mh, -eps[~part], w[~part])
        mom_dev = max(np.max(np.abs(mp["mus"] - exp_p)),
                      np.max(np.abs(mh["mus"] - exp_h)))
        A_free = profile(mp, es, exp_p) + profile(mh, -es, exp_h)
        A_free_lor = np.sum(w[None, :]*(delta/np.pi)
                            / ((es[:, None] - eps[None, :])**2 + delta**2),
                            axis=1)
        print("        U=0: max moment deviation from exact %.3e, "
              "max|A-A_exact|/max A_exact %.3e"
              % (mom_dev, np.max(np.abs(a - A_free))/np.max(A_free)),
              flush=True)

# charge gap and mu from conserved-sector ground states
gaps, mus, esec = [], [], []
for L in params["ns_gap"]:
    gl, ml, el = [], [], []
    for u in params["Us_gap"]:
        t0 = time.time()
        e = []
        for dN, dS in ((-1, -1), (0, 0), (1, 1)):
            fc = hubbard(u, shifted=False, n=L)
            fc.set_conserved_sector(Nf=L + dN, Sz=dS)
            e.append(fc.gs_energy())
        g = e[0] + e[2] - 2*e[1]
        m = 0.5*(e[2] - e[0])
        gl.append(g); ml.append(m); el.append(e)
        print("n=%d U=%.1f  gap=%.6f  mu=%.6f  (U/2=%.3f)  (%.1f s)"
              % (L, u, g, m, u/2, time.time()-t0), flush=True)
    gaps.append(gl); mus.append(ml); esec.append(el)


def lieb_wu(U, t=1.0):
    if U <= 0:
        return 0.0
    f = lambda y: np.sqrt(y*y - 1)/np.sinh(min(2*np.pi*t*y/U, 700.0))
    return 16*t*t/U*quad(f, 1, np.inf, limit=400)[0]


Ul = np.linspace(0, max(params["Us_gap"]), 200)
gap_lw = np.array([lieb_wu(u) for u in Ul])
print("free-fermion gap of the finite chain %.6f" % gap_free)
for k, u in enumerate(params["Us_gap"]):
    print("U=%.1f  finite %s  Lieb-Wu %.5f"
          % (u, " ".join("%.5f" % g[k] for g in gaps), lieb_wu(u)))

common.save(__file__, params, Us_spec=np.array(params["Us_spec"]), es=es,
            A=np.array(A_spec), sumrule=np.array(sumrule), A_free=A_free,
            A_free_lor=A_free_lor, mom_dev=np.array(mom_dev),
            Us_gap=np.array(params["Us_gap"]), ns_gap=np.array(params["ns_gap"]),
            gap=np.array(gaps),
            mu=np.array(mus), e_sectors=np.array(esec), Ul=Ul, gap_lw=gap_lw,
            gap_free=np.array(gap_free))
