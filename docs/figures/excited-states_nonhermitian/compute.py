"""Figure 5, excited-states_nonhermitian: the many-body skin effect of the
interacting Hatano-Nelson chain, from NH-DMRG.

Run with  bash docs/figures/run.sh excited-states_nonhermitian

H = sum_i [(t+g) c^dag_i c_{i+1} + (t-g) c^dag_{i+1} c_i]
    + V sum_i (n_i-1/2)(n_{i+1}-1/2),  open chain, t=1, V=1.

The density is the right-vector one, <psi_R|n_i|psi_R>/<psi_R|psi_R>.

Reference. With open boundaries H = S H0 S^-1, where H0 is the same chain
with the symmetric hopping t' = sqrt((t+g)(t-g)) and S = prod_j r^(j n_j),
r = t'/(t+g). So psi_R = S psi0, with psi0 the ground state of the
Hermitian H0, and E is the Hermitian ground-state energy. The reference
density at the full length is computed this way: an ordinary DMRG ground
state of H0, then S applied site by site (one onsite operator per site,
centred as r^(j-(n-1)/2) to keep the numbers of order one). The check of
NH-DMRG against ED is done on a shorter chain, n_ed sites.

Size. The right ground state carries the weights r^(sum_j j n_j), which
span a range of about r^(-n^2/4) across the Fock space, so the problem
gets badly conditioned fast with n and g. A first run at n=20 converged
at g=0 and 0.2, but at g=0.4 and 0.6 NH-DMRG did not reach the residual
tolerance in 5 tries (best relative residual 0.06 and 0.10, about 30
minutes each) and returned an energy below the true ground state
(-13.820 against -13.653 at g=0.4). The library warns in that case. So
the figure uses n=16, and a solve that fails the check is flagged in
`converged` and left out of the plot.

Arrays in data.npz:
    gamma      (k,)        non-reciprocity g, in units of t
    dens       (k,n)       NH-DMRG right-vector density, n sites
    dens_ref   (k,n)       reference density, Hermitian DMRG + S, n sites
    dens_bi    (k,n)       biorthogonal density <psi_L|n_i|psi_R>, n sites
    energy     (k,)        NH-DMRG energy (real part), n sites
    energy_ref (k,)        Hermitian DMRG energy of H0, n sites
    residual   (k,)        ||H psi_R - E psi_R|| / (||psi_R|| |E|), n sites
    converged  (k,)        |E_NH - E_H0| < 1e-6 and residual < 1e-4: the
                           NH-DMRG solve at n sites found the ground state
    dens_small (k,n_ed)    NH-DMRG right-vector density, n_ed sites
    dens_ed    (k,n_ed)    ED right-eigenvector density, n_ed sites
    params                 JSON string with every parameter (see common.py)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import fermionchain

params = common.start(dict(
    model="interacting Hatano-Nelson chain, spinless fermions, open",
    t=1.0,
    V=1.0,
    gamma=[0.0, 0.1, 0.2, 0.3, 0.4],
    n=16,
    n_ed=10,
    maxm=80,
    nsweeps=20,
    cutoff=1e-12,
    tol=1e-5,
    ntries=3,
    itensor_version=3,
    seeds="none (ITensor v3 random initial MPS)",
))
t, V = params["t"], params["V"]


def chain(n, tl, tr):
    """Hopping tl on c^dag_i c_{i+1}, tr on c^dag_{i+1} c_i."""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=params["itensor_version"])
    h = 0
    for i in range(n-1):
        h = h + tl*fc.Cdag[i]*fc.C[i+1] + tr*fc.Cdag[i+1]*fc.C[i]
        h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    fc.maxm = params["maxm"]
    fc.nsweeps = params["nsweeps"]
    fc.cutoff = params["cutoff"]
    return fc


def density(fc, wf, n, **kw):
    norm = fc.aMb(wf, fc.Id, wf, **kw).real
    return np.array([fc.aMb(wf, fc.N[i], wf, **kw).real/norm
                     for i in range(n)])


def nh_solve(n, g):
    fc = chain(n, t+g, t-g)
    e, psil, psir = fc.nhdmrg(tol=params["tol"], ntries=params["ntries"])
    dens = density(fc, psir, n)
    bi = np.array([(fc.aMb(psil, fc.N[i], psir)/fc.overlap(psil, psir)).real
                   for i in range(n)])
    hpsi = fc.applyoperator(fc.hamiltonian, psir)
    res = hpsi - e*psir
    resid = np.sqrt(abs(res.overlap(res))/abs(psir.overlap(psir)))/abs(e)
    return e, dens, bi, resid


def reference(n, g):
    tp = np.sqrt((t+g)*(t-g))
    r = tp/(t+g)
    fh = chain(n, tp, tp)
    e0 = fh.gs_energy()
    w = fh.get_gs()
    for j in range(n):
        w = fh.applyoperator(fh.Id + (r**(j-(n-1)/2.)-1.0)*fh.N[j], w)
    return e0, density(fh, w, n)


n, ned = params["n"], params["n_ed"]
out = {k: [] for k in ["dens", "dens_ref", "dens_bi", "energy", "energy_ref",
                       "residual", "dens_small", "dens_ed", "converged"]}
for g in params["gamma"]:
    t0 = time.time()
    # the check against ED, n_ed sites
    e_s, d_s, _, r_s = nh_solve(ned, g)
    fe = chain(ned, t+g, t-g)
    es, ws = fe.get_excited_states(mode="ED", n=1)
    d_ed = density(fe, ws[0], ned, mode="ED")
    print("g=%.2f n=%d  E_NH=%.10f  E_ED=%.10f  max|n_NH-n_ED|=%.2e  resid=%.1e"
          % (g, ned, e_s.real, np.real(es[0]), np.max(np.abs(d_s-d_ed)), r_s),
          flush=True)
    # the full length, against the similarity transform
    e, d, bi, resid = nh_solve(n, g)
    e0, dref = reference(n, g)
    print("g=%.2f n=%d  E_NH=%.10f%+.1ei  E_H0=%.10f  max|n_NH-n_ref|=%.2e  "
          "resid=%.1e  max|n_bi-1/2|=%.1e  (%.1f s)"
          % (g, n, e.real, e.imag, e0, np.max(np.abs(d-dref)), resid,
             np.max(np.abs(bi-0.5)), time.time()-t0), flush=True)
    print("   density:", np.round(d, 4), flush=True)
    conv = abs(e.real - e0) < 1e-6 and resid < 1e-4
    print("   converged:", conv, flush=True)
    for k, v in zip(out, [d, dref, bi, e.real, e0, resid, d_s, d_ed, conv]):
        out[k].append(v)

common.save(__file__, params, gamma=np.array(params["gamma"]),
            **{k: np.array(v) for k, v in out.items()})
