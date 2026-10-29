"""Screening of an impurity spin by the bath, from the ground state alone
(the NRG-MPS diagnostics of Stocker and Zilberberg, arXiv:2501.10746).

With P^sigma = n_sigma (1 - n_-sigma) the projector on the singly occupied
impurity with spin sigma, the paper defines

    M_imp = sum_sigma <P^s Sz_imp P^s>/<s|Sz|s>
    M_n   = sum_sigma <P^s Sz_n P^s>/<s|Sz|s>,   M_env = sum_n M_n

Since P^s commutes with Sz_n and P^up - P^dn = 2 Sz_imp, these are spin
correlators: M_imp = 4<Sz_imp^2>, the probability of single occupancy, and
M_n = 4<Sz_imp Sz_n>. A screened moment has M_imp close to 1 and
M_env = -M_imp. Note that M_env = -M_imp holds exactly in any ground state
with total Sz = 0, so in zero field the information is in M_imp, in how
M_n spreads over the chain and in how both respond to a field. In a field
M_n is not a connected correlator: it contains <Sz_imp><Sz_n>, so once the
field polarizes the bath (anderson_chain shifts the bath bands as well)
M_env turns positive instead of going to zero.

The Schmidt coefficients of the ground state projected on a singly
occupied impurity are the square roots of the eigenvalues of the normalized
single-occupancy block of the impurity density matrix. For one level it is
2x2, rho_ss' = <c^dag_s' c_s> with diagonal <P^s>: two equal values
1/sqrt(2) for a Kondo singlet, one value 1 for a polarized moment. For L
levels it is 2L x 2L, over the states with exactly one electron on the
impurity, rho_ab = <c^dag_b c_a prod_{c != a,b} (1 - n_c)>, and an SU(2L)
Kondo singlet gives 2L equal values 1/sqrt(2L).

For several levels, M_imp and M_n are given per level l, with Sz_l in
place of Sz_imp, and M_env per level and per channel when the bath is
given as a list of channels.

All of them are expectation values, so they work on every backend and
with mode="ED" for a cross-check. On DMRG, a bond dimension that is too
small gives a state that breaks spin symmetry; in zero field this shows up
as M_env != -M_imp and as unequal Schmidt coefficients, the same signature
as an unscreened moment, so check both against a larger maxm.
"""

import numpy as np


def kondo_screening(self, imp=None, bath=None, mode="DMRG"):
    """Return a dict with M_imp, M_env, M_n, schmidt (descending) and rho
    (the normalized single-occupancy block, modes ordered level by level,
    up then down).

    imp is an impurity orbital or a list of them (one per level), bath a
    list of bath orbitals or a list of channels, each a list of orbitals.
    Both default to the layout anderson_chain stores (fc.impurity and
    fc.channels), and otherwise to orbital 0 and all the other orbitals.
    For a single orbital M_imp and M_env are numbers and M_n an array over
    the bath; for a list, M_imp[l], M_n[l] and M_env[l] (M_env[l][k] per
    channel k when bath is a list of channels)."""
    if imp is None:
        imp = getattr(self, "impurity", [0])
        imp = imp[0] if len(imp) == 1 else imp
    single = np.ndim(imp) == 0
    imps = list(np.atleast_1d(imp))
    if bath is None:
        bath = getattr(self, "channels",
                       [[i for i in range(len(self.Nup)) if i not in imps]])
        bath = bath[0] if len(bath) == 1 else bath
    channels = bath if np.ndim(bath[0]) > 0 else None
    sites = [i for c in bath for i in c] if channels is not None else list(bath)
    vev = lambda op: self.vev(op, mode=mode)
    mimp, mn = [], []
    for i in imps:
        pup = self.Nup[i] - self.Nup[i]*self.Ndn[i] # P^up
        pdn = self.Ndn[i] - self.Nup[i]*self.Ndn[i] # P^dn
        mimp.append(np.real(vev(pup) + vev(pdn)))
        mn.append(np.array([2.0*np.real(vev((pup - pdn)*self.Sz[j]))
                            for j in sites]))
    if channels is None:
        menv = [np.sum(m) for m in mn]
    else:
        cut = np.cumsum([0] + [len(c) for c in channels])
        menv = [np.array([np.sum(m[cut[k]:cut[k+1]])
                for k in range(len(channels))]) for m in mn]
    rho = single_occupancy_block(self, imps, mode=mode)
    w = np.real(np.trace(rho))
    if w > 0.0:
        rho = rho/w
        lam = np.linalg.eigvalsh(0.5*(rho + rho.conj().T))[::-1]
        schmidt = np.sqrt(np.clip(lam, 0.0, None))
    else: # no single occupancy, nothing to decompose
        schmidt = np.zeros(len(rho))
    if single:
        mimp, menv, mn = mimp[0], menv[0], mn[0]
    return {"M_imp": mimp, "M_env": menv, "M_n": mn,
            "schmidt": schmidt, "rho": rho}


def single_occupancy_block(self, imps, mode="DMRG"):
    """<a|rho_imp|b> on the states with one electron on the impurity
    orbitals imps, modes ordered level by level, up then down"""
    c = [op for i in imps for op in (self.Cup[i], self.Cdn[i])]
    cd = [op for i in imps for op in (self.Cdagup[i], self.Cdagdn[i])]
    nm = [op for i in imps for op in (self.Nup[i], self.Ndn[i])]
    m = len(c)
    rho = np.zeros((m, m), dtype=complex)
    for a in range(m):
        for b in range(m):
            op = cd[b]*c[a] if a != b else nm[a]
            for k in range(m):
                if k != a and k != b:
                    op = op - op*nm[k] # times (1 - n_k)
            rho[a, b] = self.vev(op, mode=mode)
    return rho


class KondoScreeningMixin():
    """get_kondo_screening for a chain with spinful orbitals"""
    def get_kondo_screening(self, imp=None, bath=None, mode="DMRG"):
        """Screening diagnostics of the impurity orbital imp in the ground
        state: M_imp, M_env, M_n, Schmidt coefficients. See
        kondoscreening.py"""
        return kondo_screening(self, imp=imp, bath=bath, mode=mode)
