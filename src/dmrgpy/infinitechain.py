"""Infinite_Many_Body_Chain: a translationally-invariant chain defined by a
single n_uc-site unit cell (rather than a fixed total length), solved via
infinite DMRG (iDMRG, see pyitensor/idmrg.py for the algorithm). This is a
deliberately independent object, NOT a Many_Body_Chain subclass --
Many_Body_Chain's whole design (mode dispatch between DMRG/ED, a fixed
`self.ns`-length site list, KPM/dynamics/entanglement/excited-states
machinery, ...) assumes a *finite* chain throughout, and none of the finite
notions (an ED cross-check, a fixed number of sites) have any meaning for
an infinite system. Only `itensor_version="python"` is supported (iDMRG has
no C++/Julia backend yet) -- passing anything else raises NotImplementedError
rather than silently doing something else.

== Hamiltonian specification: L/C/R-suffixed operators ==

Every named operator is exposed three times, as flat suffixed attributes:
`SxC[i]` (site i of the *central* cell, i=0..n_uc-1 -- ordinary intra-cell
use, matching today's `Spin_Chain.Sx[i]`-style convention), `SxR[i]` (site i
of the *next* cell), and `SxL[i]` (site i of the *previous* cell, provided
purely as an ergonomic mirror so a coupling can be phrased in whichever
direction is physically natural, e.g. `SxL[i]*SxC[j]` instead of manually
reflecting indices into a `SxC[i]*SxR[j]` term).

A term is valid iff the set of cell-groups (L/C/R) it touches is a subset of
{L,C} or of {C,R} -- i.e. it may touch one or two *adjacent* cells, but
never L and R at once (that would span three cells, out of scope) --
`set_hamiltonian` validates this and raises ValueError otherwise. Any term
touching an L-suffixed operator is canonicalized by shifting it one cell to
the right (L->C, C->R), so the stored Hamiltonian is uniformly "intra-cell C
terms" + "C-to-R inter-cell terms", each physical bond attributed to
exactly one cell -- avoiding any double-counting between a cell's own C-R
terms and its right neighbor's L-C terms, the same way a user is already
expected not to redundantly write both `Sx[i]*Sx[i+1]` and `Sx[i+1]*Sx[i]`
for one bond in existing finite-chain code today.
"""

from . import multioperator
from .multioperator import MultiOperator


def _term_groups(term, n_uc):
    """The set of cell-groups ('L','C','R') a term's site indices touch,
    using this module's convention: site<0 -> L, 0<=site<n_uc -> C,
    site>=n_uc -> R."""
    groups = set()
    for _name, site in term[1:]:
        if site < 0:
            groups.add("L")
        elif site < n_uc:
            groups.add("C")
        else:
            groups.add("R")
    return groups


def _wrap_terms(op_terms, name):
    out = MultiOperator.__new__(MultiOperator)
    out.name = name
    out.op = op_terms
    out.i = len(op_terms) - 1
    return out


def _canonicalize_hamiltonian(h, n_uc):
    """Validate and split a user-supplied unit-cell MultiOperator (built
    from SxL/SxC/SxR-style operators) into (h_intra, h_inter): h_intra
    holds every term touching only the central cell (site indices
    0..n_uc-1), h_inter holds every term touching the next cell too (some
    site index >= n_uc, and at least one site index < n_uc -- see below)
    -- see this module's docstring for why L is canonicalized away
    (shifted into C/R) before this split, and why that correctly avoids
    double-counting any bond.

    A term touching *only* R-range sites (no C touch at all -- e.g. a bare
    `SxR[i]*SxR[j]` coupling) is a pure intra-cell term of the *next* cell,
    written via R purely for the user's own convenience; it is shifted
    down by n_uc and filed under h_intra too, so idmrg.py's automaton
    builder (which requires every 2-site term's smaller site index to be
    < n_uc) never has to special-case it."""
    intra_terms, inter_terms = [], []
    for term in h.op:
        distinct_sites = {site for _name, site in term[1:]}
        if len(distinct_sites) > 2:
            raise ValueError(
                "Infinite_Many_Body_Chain.set_hamiltonian: a term touches "
                "more than 2 distinct sites ({}) -- only 1- and 2-site "
                "terms are supported".format(term))
        groups = _term_groups(term, n_uc)
        if "L" in groups and "R" in groups:
            raise ValueError(
                "Infinite_Many_Body_Chain.set_hamiltonian: a term spans "
                "both the previous (L) and next (R) unit cell at once "
                "({}) -- only couplings between at most two adjacent unit "
                "cells are supported".format(term))
        if "L" in groups:
            # groups is a subset of {L,C} here (L+R together already
            # raised above) -- shifting every site by +n_uc moves an
            # L-range site (-n_uc<=site<0) into the C range and a
            # C-range site (0<=site<n_uc) into the R range, so *both*
            # labels advance one step, not just L. Leaving an original
            # "C" labeled "C" post-shift mis-filed e.g. SxL[0]*SxC[0]
            # (n_uc=2) under h_intra as ['Sx',0],['Sx',2] -- a term that
            # actually touches the next cell (site 2 >= n_uc) -- instead
            # of h_inter, violating this function's own documented split.
            term = [term[0]] + [[name, site + n_uc] for name, site in term[1:]]
            groups = {{"L": "C", "C": "R"}[g] for g in groups}
        if groups == {"R"}:
            term = [term[0]] + [[name, site - n_uc] for name, site in term[1:]]
            groups = {"C"}
        if groups <= {"C"}:
            intra_terms.append(term)
        else:
            inter_terms.append(term)
    return _wrap_terms(intra_terms, "h_intra"), _wrap_terms(inter_terms, "h_inter")


class Infinite_Many_Body_Chain:
    """iDMRG-only counterpart of Many_Body_Chain (manybodychain.py) -- see
    this module's docstring. Exposes only the narrow subset of modes v1
    supports: `set_hamiltonian`, `gs_energy` (converged energy *per site*),
    and static one-/two-point expectation values (`vev`/`correlator`)."""

    def __init__(self, site_types, itensor_version="python"):
        if itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain: only itensor_version=\"python\" "
                "is implemented (iDMRG has no C++/Julia backend yet)")
        self.itensor_version = itensor_version
        self.site_types = list(site_types)
        self.n_uc = len(self.site_types)
        if self.n_uc < 1:
            raise ValueError("Infinite_Many_Body_Chain: the unit cell needs "
                              "at least one site")
        if self.n_uc > 2:
            raise NotImplementedError(
                "Infinite_Many_Body_Chain: unit cells with more than 2 "
                "sites (n_uc={}) are not supported yet -- the growth loop's "
                "micro-step sublattice pairing (pyitensor/idmrg.py's "
                "p_L=mstep, p_R=n_uc-1-mstep) only produces genuinely "
                "adjacent active sites for n_uc in {{1, 2}}; for n_uc>=3 it "
                "pairs non-adjacent sublattice positions whose automaton "
                "channels carry different, bond-specific pending content "
                "(not just a matching channel count), so contracting them "
                "together would silently mix unrelated bond content rather "
                "than just fail to contract. Fixing this needs a genuine "
                "redesign of the micro-step growth scheme, not a small "
                "patch -- flagged as a known follow-up rather than "
                "attempted here.".format(self.n_uc))

        self.maxm = 30          # wavefunction bond dimension cap
        self.cutoff = 1e-12     # SVD truncation
        self.maxiter = 200      # iDMRG macro-iterations (growth steps)
        self.etol = 1e-10       # energy-density convergence tolerance
        self.niter = 30         # per-micro-step Lanczos iteration count --
                                 # note idmrg.py's _local_two_site_solve
                                 # always runs at least 200 (a validated,
                                 # load-bearing floor, see its own comment),
                                 # so any value <=200 here has no effect;
                                 # only raising it above 200 does anything
        self.verbose = False

        self.hamiltonian = None  # user-facing MultiOperator (L/C/R indices)
        self._h_intra = None
        self._h_inter = None
        self._result = None      # pyitensor.idmrg.IDMRGResult once converged
        self.e0 = None           # converged ground-state energy per site
        self.converged = None

    def get_operator(self, name, i, group="C"):
        """A bare, symbolic 1-site MultiOperator for `name` at site `i`
        (0..n_uc-1) of cell-group `group` ('L','C','R') -- mirrors
        Many_Body_Chain.get_operator's `multioperator.obj2MO([[name,i]])`,
        just with the L/C/R site-index offset applied first."""
        if group == "C":
            site = i
        elif group == "R":
            site = self.n_uc + i
        elif group == "L":
            site = i - self.n_uc
        else:
            raise ValueError("get_operator: group must be 'L', 'C' or 'R', got {!r}".format(group))
        return multioperator.obj2MO([[name, site]])

    def set_hamiltonian(self, h):
        """h: a MultiOperator built from this chain's SxL/SxC/SxR-style
        (or get_operator(...,group=...)-built) operators -- see this
        module's docstring for the validity/canonicalization rules."""
        # NOTE: v1 is Hermitian-only (see this module's docstring) and
        # idmrg_ground_state has no Hermiticity check of its own, unlike
        # the finite-chain path (groundstate.py), which gates on
        # MultiOperator.is_hermitian() and routes a non-Hermitian
        # Hamiltonian to a dedicated NH-DMRG solver instead -- a real gap
        # for this module. A symbolic is_hermitian() check was tried here
        # and reverted: confirmed directly (also on an ordinary *finite*
        # Spin_Chain, so this is a pre-existing, general limitation, not
        # specific to this module) that is_hermitian()'s simplify() step
        # does not recognize that reversing a product of operators on
        # *different* sites (as get_dagger() does) gives an equivalent
        # term to the un-reversed original, so it false-rejects an
        # ordinary Sx[i]*Sx[j]+Sy[i]*Sy[j]+Sz[i]*Sz[j] Heisenberg term --
        # the single most common Hamiltonian shape this module exists to
        # support. Left as a known, documented gap rather than shipping
        # that regression.
        h = multioperator.obj2MO(h)
        self._h_intra, self._h_inter = _canonicalize_hamiltonian(h, self.n_uc)
        self.hamiltonian = h
        self._result = None
        self.e0 = None
        self.converged = None

    def gs_energy(self):
        """Run iDMRG to convergence (or self.maxiter macro-iterations) and
        return the ground-state energy *per site* -- the physically
        meaningful observable for an infinite system (a total energy would
        be unboundedly large)."""
        if self._h_intra is None:
            raise RuntimeError(
                "Infinite_Many_Body_Chain.gs_energy called before set_hamiltonian")
        from .pyitensor import idmrg
        self._result = idmrg.idmrg_ground_state(
            self.site_types, self._h_intra.op, self._h_inter.op, self.n_uc,
            maxm=self.maxm, cutoff=self.cutoff, maxiter=self.maxiter,
            etol=self.etol, niter=self.niter, verbose=self.verbose)
        self.e0 = self._result.e0
        self.converged = self._result.converged
        return self.e0

    def vev(self, opname, p, group="C"):
        """<opname> at site p (0..n_uc-1) of cell-group `group` -- by
        translational invariance this is identical for every group, `group`
        is accepted purely so callers can write whichever reads most
        naturally (SxL/SxC/SxR all describe the same infinite chain)."""
        if group not in ("L", "C", "R"):
            raise ValueError("vev: group must be 'L', 'C' or 'R', got {!r}".format(group))
        if not (0 <= p < self.n_uc):
            raise ValueError("vev: p must be in 0..{} (n_uc-1), got {!r}".format(
                self.n_uc - 1, p))
        if self._result is None:
            self.gs_energy()
        from .pyitensor import idmrg
        return idmrg.onsite_expectation(self._result, opname, p)

    def correlator(self, opname_i, p_i, opname_j, r):
        """<opname_i(site p_i) opname_j(site p_i + r)>, r measured in
        physical sites (r>=0) -- see pyitensor.idmrg.two_point_correlator
        for the r=0 (same-site) convention."""
        if not (0 <= p_i < self.n_uc):
            raise ValueError("correlator: p_i must be in 0..{} (n_uc-1), got {!r}".format(
                self.n_uc - 1, p_i))
        if self._result is None:
            self.gs_energy()
        from .pyitensor import idmrg
        return idmrg.two_point_correlator(self._result, opname_i, p_i, opname_j, r)


class Infinite_Spin_Chain(Infinite_Many_Body_Chain):
    """Infinite counterpart of spinchain.Spin_Chain: a unit cell of `spins`
    (a list of n_uc spin labels, e.g. "1/2","1",... -- see
    spinchain.label2site), with SxC/SyC/SzC/SxL/SyL/SzL/SxR/SyR/SzR
    per-site operator lists (length n_uc each)."""

    def __init__(self, spins, **kwargs):
        from .spinchain import label2site
        site_types = [label2site[s] for s in spins]
        Infinite_Many_Body_Chain.__init__(self, site_types, **kwargs)
        n = self.n_uc
        for name in ("Sx", "Sy", "Sz"):
            for group in ("L", "C", "R"):
                setattr(self, name + group,
                        [self.get_operator(name, i, group) for i in range(n)])
