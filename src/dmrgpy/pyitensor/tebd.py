"""TEBD: 2nd-order Trotter (odd/even bond) real-time evolution for
strictly nearest-neighbor Hamiltonians.

Unlike tdvp.py's per-bond, per-step Krylov exponentiation of an
environment-projected *effective* Hamiltonian, each bond's evolution gate
exp(-i*tau*h_bond) here is the exact exponential of the *bare* local
2-site Hamiltonian (a small dense matrix, exponentiated once via
scipy.linalg.expm) -- the standard Vidal TEBD algorithm. Since H is
time-independent for every caller in this codebase (a "quench" holds H
fixed across all nt steps), every gate is built once at setup and reused
unchanged for every subsequent time step: no per-step diagonalization at
all, only a gate contraction + SVD truncation per bond per step. This is
what makes TEBD cheaper than TDVP for models where it applies, at the
cost of only supporting strictly nearest-neighbor H (any term spanning 3+
distinct sites raises NotImplementedError -- long-range Hamiltonians need
TDVP or the MPO-Taylor path instead).

One-site (onsite) terms are split half-and-half onto each of a site's
neighboring bonds (full weight at a chain boundary, where a site has only
one neighbor) -- the same convention TeNPy's NearestNeighborModel uses.
Bare scalar terms (no site reference at all, e.g. a ground-state-energy
shift folded in as `coef*Identity`) are absorbed entirely into bond 1;
since they act as Identity on every other bond too, which bond carries
them doesn't matter.

Fermionic terms need no special-casing here beyond what HTerm.resolve()
already does: for any *physical* (fermion-parity-conserving) Hamiltonian
term restricted to two adjacent sites, resolve()'s Jordan-Wigner carry
always returns to its initial (False) parity by the far site, so every
site beyond the term's own span gets a plain Id, not a trailing F string --
i.e. embedding kron(mats[lo-1], mats[hi-1]) with implicit identity
elsewhere reproduces exactly the same operator the true (JW-threaded) MPO
Hamiltonian represents. bond_hamiltonians() asserts this explicitly rather
than trusting it silently, so a non-parity-conserving term (which would
need a JW string extending past the two-site bond, incompatible with a
local gate) fails loudly instead of silently producing wrong physics.
"""

import numpy as np
from scipy.linalg import expm

from .autompo import AutoMPO
from .mpscontainer import _link_at
from .svd import svd
from .tensor import ITensor, noPrime


def _true_span(mats, dims):
    """The 1-based (lo, hi) site range where `mats` (one full-chain matrix
    per site, from HTerm.resolve()) differs from the local identity --
    i.e. the sites this term actually, non-trivially acts on.

    Determining this from the *resolved* matrices rather than the raw
    op-name list is necessary, not just more convenient: two independent
    sources produce op-name entries that look like they reach a far-away
    site but algebraically compose to a no-op there. (1)
    MultiOperator.identity() (multioperator.py) hardcodes its bare "Id"
    factor at a fixed placeholder site regardless of context, so e.g.
    `N[i]-0.5` (a common particle-hole-symmetric-offset pattern used
    throughout dmrgpy's own Hubbard/UV-chain examples) expands to terms
    like [('N',i),('Id',<placeholder>)]. (2) jordanwigner_spinful.py's
    explicit Jordan-Wigner threading for native spinful (Hubbard/Electron)
    sites can insert paired F-string factors at a site that multiply out
    to the identity (F @ F = Id) rather than cancelling them away
    algebraically first -- confirmed directly: a plain nearest-neighbor
    hopping term on a Spinful_Fermionic_Chain_Native chain carries extra
    ('F', site) factors reaching back to site 1 in its raw op list, even
    though the resolved matrix at site 1 is exactly the identity. Treating
    either of these as "real" span in a naive op-name scan would flag
    strictly local physics as spurious long-range terms. This also
    subsumes the previous separate fermion-parity-conservation check
    (an unresolved JW carry manifests as a non-identity matrix at every
    site through the end of the chain, so it is caught by the same
    hi-lo>1 test below, just without a special-cased second pass).
    """
    touched = [i for i, m in enumerate(mats) if not np.allclose(m, np.eye(dims[i], dtype=complex))]
    if not touched:
        return None
    return touched[0] + 1, touched[-1] + 1


def bond_hamiltonians(sites, terms):
    """terms: dmrgpy MultiOperator.to_terms() shape, i.e.
    [(coef, [(opname,site), ...]), ...]. Returns {b: (D,D) ndarray} for
    b=1..n-1 (the bond between sites b and b+1), in the same standard
    (out,in) convention as HTerm.resolve()."""
    n = sites.length()
    if n < 2:
        raise ValueError("bond_hamiltonians: TEBD needs at least 2 sites")
    ampo = AutoMPO.from_terms(sites, terms)
    dims = [sites.dim(i) for i in range(1, n + 1)]
    H = {b: np.zeros((dims[b - 1] * dims[b], dims[b - 1] * dims[b]), dtype=complex)
         for b in range(1, n)}

    for term in ampo.terms:
        mats = term.resolve(sites)  # standard (out,in) per-site matrices
        span = _true_span(mats, dims)
        if span is None:
            H[1] += term.coef * np.eye(dims[0] * dims[1], dtype=complex)
            continue
        lo, hi = span
        if hi - lo > 1:
            raise NotImplementedError(
                "TEBD requires a nearest-neighbor Hamiltonian; found a term "
                "with non-trivial support spanning sites {}..{}".format(lo, hi))
        if hi == lo:
            local = term.coef * mats[lo - 1]
            targets = []
            if lo > 1:
                targets.append((lo - 1, "right"))  # bond (lo-1,lo): term is the right site
            if lo < n:
                targets.append((lo, "left"))       # bond (lo,lo+1): term is the left site
            weight = 0.5 if len(targets) == 2 else 1.0
            for b, side in targets:
                if side == "right":
                    other_dim = dims[b - 1]
                    piece = np.kron(np.eye(other_dim, dtype=complex), local)
                else:
                    other_dim = dims[b]
                    piece = np.kron(local, np.eye(other_dim, dtype=complex))
                H[b] += weight * piece
        else:
            H[lo] += term.coef * np.kron(mats[lo - 1], mats[hi - 1])
    return H


def _bond_gate(sites, i, h_mat, tau):
    """The 2-site ITensor gate exp(-i*tau*h_mat) for bond (i,i+1)."""
    s_i, s_j = sites.si(i), sites.si(i + 1)
    d_i, d_j = s_i.dim, s_j.dim
    U = expm(-1j * tau * h_mat)
    stored = U.T.reshape(d_i, d_j, d_i, d_j)  # (in_i,in_j,out_i,out_j)
    return ITensor((s_i, s_j, s_i.prime(1), s_j.prime(1)), stored)


def _apply_bond_gate(psi, i, gate, cutoff, maxdim):
    """Contract `gate` onto bond (i,i+1) and SVD-truncate. Canonicalizes
    psi's orthogonality center to site i first (via position(), the same
    primitive autompo.py/mpsalgebra.py use for their own truncating
    sweeps), which is what makes the local SVD truncation-optimal
    regardless of which bond was last touched. Mutates psi in place."""
    psi.position(i, cutoff=cutoff, maxdim=maxdim)
    left_link = _link_at(psi, i, i - 1)
    s_i = next(ind for ind in psi.A(i).inds if ind.hastags("Site"))
    theta = noPrime(gate * (psi.A(i) * psi.A(i + 1)), "Site")
    left_inds = ([left_link] if left_link else []) + [s_i]
    U, S, V, _spec = svd(theta, left_inds, cutoff=cutoff, maxdim=maxdim)
    psi.set_A(i, U)
    psi.set_A(i + 1, S * V)
    psi.center = i + 1


class TEBDEvolver:
    """Precomputes every bond's dt/2 and dt evolution gates once from a
    (static, nearest-neighbor) Hamiltonian, then applies a 2nd-order
    Trotter step exp(-i dt/2 H_odd) exp(-i dt H_even) exp(-i dt/2 H_odd)
    per call to step() -- H_odd/H_even being the sums over odd-indexed
    and even-indexed bonds respectively, each internally commuting since
    every bond in one group acts on disjoint sites."""

    def __init__(self, sites, terms, dt, cutoff, maxdim):
        self.sites = sites
        self.n = sites.length()
        self.cutoff = cutoff
        self.maxdim = maxdim
        h_bonds = bond_hamiltonians(sites, terms)
        self._odd = list(range(1, self.n, 2))
        self._even = list(range(2, self.n, 2))
        # step() only ever applies the half-dt gate on odd bonds and the
        # full-dt gate on even bonds -- build only those, not both gates
        # for every bond.
        self._gates_half = {b: _bond_gate(sites, b, h_bonds[b], dt / 2.0) for b in self._odd}
        self._gates_full = {b: _bond_gate(sites, b, h_bonds[b], dt) for b in self._even}

    def step(self, psi):
        """One full time step of size dt. Mutates psi in place, returns it."""
        for b in self._odd:
            _apply_bond_gate(psi, b, self._gates_half[b], self.cutoff, self.maxdim)
        for b in self._even:
            _apply_bond_gate(psi, b, self._gates_full[b], self.cutoff, self.maxdim)
        for b in self._odd:
            _apply_bond_gate(psi, b, self._gates_half[b], self.cutoff, self.maxdim)
        return psi
