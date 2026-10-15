# Mott gap of the half-filled Hubbard chain. Local spectral function A_ii(w)
# of a bulk site (spin up) from KPM, and the charge gap from the ground
# states of the sectors N-1, N, N+1, against the Lieb-Wu gap.
import numpy as np
from scipy.integrate import quad
from dmrgpy import fermionchain, kpmdmrg

n = 12
site = n//2
Us_spec = [0.0, 2.0, 4.0, 8.0]
Us_gap = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0]
ns_gap = [n, 2*n]
delta = 0.25
es = np.linspace(-8.0, 8.0, 800)


def hubbard(U, n, shifted):
    # shifted: U (n_up-1/2)(n_dn-1/2), particle-hole symmetric, so mu=0
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = sum(-cd[i]*c[i+1] - cd[i+1]*c[i] for i in range(n-1)
            for cd, c in ((fc.Cdagup, fc.Cup), (fc.Cdagdn, fc.Cdn)))
    if shifted:
        h = h + sum(U*(fc.Nup[i] - 0.5)*(fc.Ndn[i] - 0.5) for i in range(n))
    else:
        h = h + sum(U*fc.Nup[i]*fc.Ndn[i] for i in range(n))
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps, fc.cutoff, fc.kpmmaxm = 80, 12, 1e-12, 80
    return fc


def kpm(fc, A, B):
    # KPM moments of (A,B) with the rescaling needed to rebuild a spectrum
    mus, emin, emax, scale, npol, d = fc.get_dynamical_correlator_moments(
        name=(A, B), delta=delta)
    return dict(mus=np.asarray(mus), emin=emin, emax=emax, scale=scale,
                npol=npol, d=d, e0=fc.e0)


def spectrum(x, m):
    return kpmdmrg.dynamical_correlator_from_moments(
        m["mus"], m["emin"], m["emax"], m["scale"], m["npol"], x,
        kernel="jackson", delta=m["d"], origin=m["e0"])[1].real


def exact(m, poles, weights):
    # m with its moments replaced by the exact ones of a set of poles
    x = (poles + m["e0"] - (m["emin"] + m["emax"])/2)*m["scale"]
    j = np.arange(len(m["mus"]))[:, None]
    return dict(m, mus=np.sum(weights*np.cos(j*np.arccos(x)), axis=1))


# free fermions at the site: levels, weights, empty levels at half filling
eps, phi = np.linalg.eigh(-(np.eye(n, k=1) + np.eye(n, k=-1)))
w, empty = phi[site]**2, eps > 0

A = []
for U in Us_spec:
    fc = hubbard(U, n, shifted=True)
    fc.gs_energy()
    mp = kpm(fc, fc.Cup[site], fc.Cdagup[site])  # addition, w > 0
    mh = kpm(fc, fc.Cdagup[site], fc.Cup[site])  # removal, w < 0
    A.append(spectrum(es, mp) + spectrum(-es, mh))
    if U == 0.0:  # exact free-fermion moments, same rescaling and kernel
        A_free = (spectrum(es, exact(mp, eps[empty], w[empty]))
                  + spectrum(-es, exact(mh, -eps[~empty], w[~empty])))

# charge gap and mu from the sectors N-1, N, N+1, with U n_up n_dn
gap, mu = np.zeros((2, len(ns_gap), len(Us_gap)))
for a, L in enumerate(ns_gap):
    for b, U in enumerate(Us_gap):
        e = []
        for dN in (-1, 0, 1):
            fc = hubbard(U, L, shifted=False)
            fc.set_conserved_sector(Nf=L + dN, Sz=dN)
            e.append(fc.gs_energy())
        gap[a, b], mu[a, b] = e[0] + e[2] - 2*e[1], (e[2] - e[0])/2


def lieb_wu(U):
    # charge gap of the infinite chain, t=1
    f = lambda y: np.sqrt(y*y - 1)/np.sinh(min(2*np.pi*y/U, 700.0))
    return 16/U*quad(f, 1, np.inf, limit=400)[0] if U > 0 else 0.0


Ul = np.linspace(0, max(Us_gap), 200)
np.savez("data.npz", Us_spec=Us_spec, es=es, A=A, A_free=A_free,
         Us_gap=Us_gap, ns_gap=ns_gap, gap=gap, mu=mu,
         Ul=Ul, gap_lw=[lieb_wu(U) for U in Ul])
