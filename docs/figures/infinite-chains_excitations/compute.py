"""Figure 21, infinite-chains_excitations: excitation bands of infinite
chains from the quasiparticle ansatz on top of VUMPS, drawn over the
dynamical structure factor S(k,omega) of the same ansatz.

(a) transverse-field Ising chain H = -sum sx sx - g sum sz, g=2.5,
    S^xx(k,omega); exact band 2 sqrt(1+g^2-2g cos k).
(b) spin-1 Heisenberg chain, S^zz(k,omega); Haldane gap 0.4104789 at k=pi.

Run with  bash docs/figures/run.sh infinite-chains_excitations

Arrays in data.npz:
    ks_a      (na,)     momenta of (a), radians per site
    band_a    (na,)     lowest excitation energy E(k) of (a), from
                        excitation_energies(k, n=1)
    exact_a   (na,)     exact free-fermion band of (a)
    es_a      (ma,)     energies of the map of (a)
    S_a       (na, ma)  S^xx(k, omega) of (a), Lorentzian broadening delta
    e0_a      ()        VUMPS energy per site of (a)
    e0exact_a ()        exact energy per site of (a)
    conv_a    ()        ic.converged of (a)
    ks_b      (nb,)     momenta of (b)
    e_b       (nb, nbr) branch energies of (b), lowest first
    w_b       (nb, nbr) spectral weights |<k,a|S^z_k|GS>|^2 of (b)
    tot_b     (nb,)     total weight (static structure factor S^zz(k)) of (b)
    es_b      (mb,)     energies of the map of (b)
    S_b       (nb, mb)  S^zz(k, omega) of (b), broadened with delta
    e0_b      ()        VUMPS energy per site of (b)
    conv_b    ()        ic.converged of (b)
    gap_ref   ()        Haldane gap 0.4104789 (reference)
    e0ref_b   ()        energy per site -1.401484039 (reference)
    params              JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import infinitechain

params = common.start(dict(
    itensor_version="python",
    gs_method="vumps",
    a_model="TFIM, H=-4 SxSx - 2 g Sz (Pauli -sx sx - g sz), n_uc=1",
    a_g=2.5, a_maxm=4, a_maxiter=800, a_nrestarts=6,
    a_nk=int(os.environ.get("FIG_NKA", 61)),
    a_energies=[1.5, 8.5, 400], a_delta=0.08,
    b_model="spin-1 Heisenberg chain, J=1, n_uc=1",
    b_maxm=int(os.environ.get("FIG_MAXMB", 16)), b_maxiter=800,
    b_nrestarts=6, b_nk=int(os.environ.get("FIG_NKB", 97)), b_nbranch=12,
    b_energies=[0.0, 6.0, 400], b_delta=0.08,
    seeds="none (VUMPS random starts, recorded through e0 and converged)",
))
g = params["a_g"]

# ---------------- (a) transverse-field Ising chain ------------------------
t0 = time.time()
ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version="python")
ic.set_hamiltonian(-4.0*ic.SxC[0]*ic.SxR[0] - 2.0*g*ic.SzC[0])
ic.gs_method = "vumps"
ic.maxm = params["a_maxm"]
ic.maxiter = params["a_maxiter"]
ic.vumps_nrestarts = params["a_nrestarts"]
e0_a = ic.gs_energy()
conv_a = bool(ic.converged)


def eps(k):
    return 2.0*np.sqrt(1.0 + g**2 - 2.0*g*np.cos(k))


# exact ground-state energy per site: -(1/2pi) int eps(k)/2 dk
kq = np.linspace(-np.pi, np.pi, 20001)
e0exact_a = -np.mean(eps(kq[:-1]))/2   # -(1/2) <eps>, periodic mean
print("(a) e0=%.12f exact=%.12f diff=%.2e converged=%s"
      % (e0_a, e0exact_a, e0_a - e0exact_a, conv_a), flush=True)
ks_a = np.linspace(-np.pi, np.pi, params["a_nk"])
band_a = np.array([ic.excitation_energies(k, n=1)[0] for k in ks_a]).real
exact_a = eps(ks_a)
print("(a) max |E(k)-exact| = %.2e" % np.max(np.abs(band_a - exact_a)))
lo, hi, ne = params["a_energies"]
_, es_a, S_a = ic.dynamical_structure_factor(
    "Sx", ks=ks_a, energies=np.linspace(lo, hi, ne),
    delta=params["a_delta"], n=1)
print("(a) done in %.1f s" % (time.time() - t0), flush=True)

# ---------------- (b) spin-1 Heisenberg chain -----------------------------
t0 = time.time()
ich = infinitechain.Infinite_Spin_Chain(["1"], itensor_version="python")
ich.set_hamiltonian(ich.SxC[0]*ich.SxR[0] + ich.SyC[0]*ich.SyR[0]
                    + ich.SzC[0]*ich.SzR[0])
ich.gs_method = "vumps"
ich.maxm = params["b_maxm"]
ich.maxiter = params["b_maxiter"]
ich.vumps_nrestarts = params["b_nrestarts"]
e0_b = ich.gs_energy()
conv_b = bool(ich.converged)
E0REF_B, GAP_REF = -1.401484038971, 0.4104789
print("(b) e0=%.10f ref=%.10f diff=%.2e converged=%s (%.1f s)"
      % (e0_b, E0REF_B, e0_b - E0REF_B, conv_b, time.time() - t0), flush=True)
ks_b = np.linspace(0.0, 2*np.pi, params["b_nk"])
e_b, w_b, tot_b = [], [], []
for k in ks_b:
    e, w, tot = ich.spectral_weights("Sz", k, n=params["b_nbranch"],
                                     return_total=True)
    e_b.append(np.real(e))
    w_b.append(np.real(w))
    tot_b.append(np.real(tot))
e_b, w_b, tot_b = np.array(e_b), np.array(w_b), np.array(tot_b)
ipi = int(np.argmin(np.abs(ks_b - np.pi)))
print("(b) gap E(pi)=%.8f ref=%.7f diff=%.2e; triplet weight at pi %.4f"
      % (e_b[ipi, 0], GAP_REF, e_b[ipi, 0] - GAP_REF,
         w_b[ipi, :3].sum()/tot_b[ipi]))
lo, hi, ne = params["b_energies"]
es_b = np.linspace(lo, hi, ne)
dl = params["b_delta"]
lor = (dl/np.pi)/((es_b[None, :, None] - e_b[:, None, :])**2 + dl**2)
S_b = np.einsum("kea,ka->ke", lor, w_b)
print("(b) done in %.1f s" % (time.time() - t0), flush=True)

common.save(__file__, params, ks_a=ks_a, band_a=band_a, exact_a=exact_a,
            es_a=np.asarray(es_a), S_a=np.asarray(S_a).real,
            e0_a=np.array(e0_a), e0exact_a=np.array(e0exact_a),
            conv_a=np.array(conv_a), ks_b=ks_b, e_b=e_b, w_b=w_b,
            tot_b=tot_b, es_b=es_b, S_b=S_b, e0_b=np.array(e0_b),
            conv_b=np.array(conv_b), gap_ref=np.array(GAP_REF),
            e0ref_b=np.array(E0REF_B))
