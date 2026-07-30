"""Infinite_Many_Body_Chain: a translationally-invariant chain defined by a
single n_uc-site unit cell (rather than a fixed total length), solved via
infinite DMRG (iDMRG, see pyitensor/idmrg.py for the algorithm, or
mpscpp3/chain_session.h's Chain::idmrg_ground_state for the ITensor v3 C++
port of the same algorithm). This is a deliberately independent object, NOT
a Many_Body_Chain subclass -- Many_Body_Chain's whole design (mode dispatch
between DMRG/ED, a fixed `self.ns`-length site list, KPM/dynamics/
entanglement/excited-states machinery, ...) assumes a *finite* chain
throughout, and none of the finite notions (an ED cross-check, a fixed
number of sites) have any meaning for an infinite system.
`itensor_version="python"` (default) or `3` are supported -- passing
anything else raises NotImplementedError rather than silently doing
something else. The v3 C++ backend computes the energy density only (no
static-correlator support yet, see Infinite_Many_Body_Chain.gs_energy's own
comment); `vev`/`correlator` still require `itensor_version="python"`
regardless of which backend `gs_energy` itself used.

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

import numpy as np

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


def _shift_terms(op_terms, shift):
    """A copy of a MultiOperator.op-format term list with every site index
    increased by `shift` -- used by _window_hamiltonian to tile h_intra/
    h_inter (already canonicalized to a single unit cell, 0-based site
    indices) across a finite window's unit-cell repeats."""
    return [[term[0]] + [[name, site + shift] for name, site in term[1:]]
            for term in op_terms]


def _window_hamiltonian(h_intra, h_inter, n_uc, n_window):
    """A finite, open-boundary MultiOperator over n_window*n_uc sites (0-
    based), built by tiling h_intra once per cell (0..n_window-1) and
    h_inter once per *adjacent pair* of cells (0..n_window-2 -- one fewer
    bond than a periodic ring, since there is no cell n_window to couple
    the last cell's own h_inter to). See
    Infinite_Many_Body_Chain.kpm_finite's own docstring for why this
    finite tiling is used (and its approximation to the true infinite
    chain) rather than an exact infinite-size construction."""
    terms = []
    for c in range(n_window):
        terms += _shift_terms(h_intra.op, c * n_uc)
    for c in range(n_window - 1):
        terms += _shift_terms(h_inter.op, c * n_uc)
    return _wrap_terms(terms, "h_window")


class Infinite_Many_Body_Chain:
    """iDMRG-only counterpart of Many_Body_Chain (manybodychain.py) -- see
    this module's docstring. Exposes only a narrow subset of modes:
    `set_hamiltonian`, `gs_energy` (converged energy *per site*), static
    one-/two-point expectation values (`vev`/`correlator`), an
    open-boundary-window approximation to dynamical correlators
    (`kpm_finite`), and an infinite-boundary-condition, real-time-TDVP
    dynamical correlator (`td_dynamical_correlator`, itensor_version=
    "python" only -- see its own docstring for how it differs from
    `kpm_finite`). No excited states, no entanglement/entropy yet.

    `itensor_version="python"` (default) or `3` -- the ITensor v3 C++
    backend (`mpscpp3/chain_session.h`'s `Chain::idmrg_ground_state`)
    computes the energy density only; it has no static-correlator support
    yet (see `gs_energy`'s own comment), so `vev`/`correlator` still
    require `itensor_version="python"` regardless of what backend
    `gs_energy` itself was run with. `itensor_version=2` and
    `"julia_live"`/`"julia"` have no iDMRG port at all and raise
    `NotImplementedError`."""

    def __init__(self, site_types, itensor_version="python"):
        if itensor_version not in ("python", 3):
            raise NotImplementedError(
                "Infinite_Many_Body_Chain: itensor_version={!r} is not "
                "implemented -- only \"python\" and 3 (the ITensor v3 C++ "
                "backend) support iDMRG".format(itensor_version))
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
        self.niter = 30         # itensor_version="python": per-micro-step
                                 # Lanczos iteration count -- note
                                 # idmrg.py's _local_two_site_solve always
                                 # runs at least 200 (a validated,
                                 # load-bearing floor, see its own comment),
                                 # so any value <=200 here has no effect;
                                 # only raising it above 200 does anything.
                                 # itensor_version=3: reused as the local
                                 # 2-site solve's Arnoldi Krylov dimension
                                 # (Chain::idmrg_ground_state's krylovdim
                                 # argument).
        self.restarts = 2       # itensor_version=3 only: Arnoldi restart
                                 # count for the local 2-site solve
                                 # (Chain::idmrg_ground_state's restarts
                                 # argument) -- no equivalent knob on the
                                 # "python" backend's Lanczos solve.
        self.verbose = False

        self.hamiltonian = None  # user-facing MultiOperator (L/C/R indices)
        self._h_intra = None
        self._h_inter = None
        self._result = None      # pyitensor.idmrg.IDMRGResult once converged
                                  # (itensor_version="python" only -- the v3
                                  # backend has no correlator machinery, see
                                  # gs_energy/vev/correlator)
        self.e0 = None           # converged ground-state energy per site
        self.converged = None
        self._excitation_env = None  # pyitensor.idmrg_excitations.
                                      # ExcitationEnvironment, lazily built by
                                      # excitation_energies/excitation_gap and
                                      # cached across repeated calls (e.g. a
                                      # dispersion/gap k-scan) -- invalidated
                                      # in set_hamiltonian exactly like
                                      # self._result already is.
        self._session3 = None    # itensor_version=3 only: the mpscpp3
                                  # Chain instance gs_energy() last ran
                                  # idmrg_ground_state on, kept alive so
                                  # td_dynamical_correlator can reuse its
                                  # private converged-environment snapshot
                                  # (Chain::idmrg_ground_state's own
                                  # env_window_boundary-equivalent state,
                                  # never exposed to Python directly -- see
                                  # chain_session.h's own comment) -- a
                                  # *fresh* Chain() has no such snapshot,
                                  # so re-running gs_energy() on a new
                                  # instance would silently lose it.

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
        self._excitation_env = None
        self._session3 = None

    def gs_energy(self):
        """Run iDMRG to convergence (or self.maxiter macro-iterations) and
        return the ground-state energy *per site* -- the physically
        meaningful observable for an infinite system (a total energy would
        be unboundedly large).

        itensor_version="python" runs pyitensor/idmrg.py's own growing
        algorithm and keeps the resulting IDMRGResult in self._result for
        vev/correlator to reuse (see those methods). itensor_version=3
        instead calls the compiled mpscpp3 backend's
        Chain::idmrg_ground_state directly (energy density only -- no
        correlator support there yet, see this module's own docstring),
        so self._result is left None on that path; vev/correlator raise
        NotImplementedError rather than silently misusing a stale/absent
        result."""
        if self._h_intra is None:
            raise RuntimeError(
                "Infinite_Many_Body_Chain.gs_energy called before set_hamiltonian")
        self._excitation_env = None  # about to rebuild self._result; see __init__'s own comment
        if self.itensor_version == "python":
            from .pyitensor import idmrg
            self._result = idmrg.idmrg_ground_state(
                self.site_types, self._h_intra.op, self._h_inter.op, self.n_uc,
                maxm=self.maxm, cutoff=self.cutoff, maxiter=self.maxiter,
                etol=self.etol, niter=self.niter, verbose=self.verbose)
            self.e0 = self._result.e0
            self.converged = self._result.converged
        else:  # itensor_version == 3
            from . import cppext
            backend = cppext.get_backend(3)
            if backend is None:
                raise RuntimeError(
                    "Infinite_Many_Body_Chain.gs_energy: itensor_version=3 "
                    "requested but the mpscpp3 (ITensor v3) extension is "
                    "not compiled -- run install.py --itensor-version=3 "
                    "first, or use itensor_version=\"python\" instead")
            chain = backend.Chain(self.site_types)
            chain.set_verbose(self.verbose)
            terms_intra = self._h_intra.to_terms()
            terms_inter = self._h_inter.to_terms()
            density, converged, _niter_done = chain.idmrg_ground_state(
                terms_intra, terms_inter, self.maxm, self.cutoff,
                self.maxiter, self.etol, self.niter, self.restarts)
            self._session3 = chain  # kept alive for td_dynamical_correlator -- see __init__'s own comment
            self._result = None
            self.e0 = density
            self.converged = converged
        return self.e0

    def vev(self, opname, p, group="C"):
        """<opname> at site p (0..n_uc-1) of cell-group `group` -- by
        translational invariance this is identical for every group, `group`
        is accepted purely so callers can write whichever reads most
        naturally (SxL/SxC/SxR all describe the same infinite chain).

        Only supported for itensor_version="python" -- the ITensor v3 C++
        backend has no static-correlator machinery yet (see gs_energy's
        own comment)."""
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.vev: only itensor_version="
                "\"python\" supports static correlators -- run a separate "
                "itensor_version=\"python\" chain for vev/correlator, or "
                "reuse pyitensor.idmrg's correlator functions directly "
                "against a \"python\"-backend IDMRGResult")
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
        for the r=0 (same-site) convention.

        Only supported for itensor_version="python" -- see vev's own
        comment."""
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.correlator: only itensor_version="
                "\"python\" supports static correlators -- run a separate "
                "itensor_version=\"python\" chain for vev/correlator, or "
                "reuse pyitensor.idmrg's correlator functions directly "
                "against a \"python\"-backend IDMRGResult")
        if not (0 <= p_i < self.n_uc):
            raise ValueError("correlator: p_i must be in 0..{} (n_uc-1), got {!r}".format(
                self.n_uc - 1, p_i))
        if self._result is None:
            self.gs_energy()
        from .pyitensor import idmrg
        return idmrg.two_point_correlator(self._result, opname_i, p_i, opname_j, r)

    def _get_excitation_environment(self):
        """Lazily build (or reuse) the pyitensor.idmrg_excitations.
        ExcitationEnvironment for the converged ground state -- expensive
        (a null-space computation plus two dense linear solves), so it is
        cached on self._excitation_env across repeated excitation_energies/
        excitation_gap calls (e.g. a k-scan), and invalidated whenever
        set_hamiltonian/gs_energy reruns (see __init__'s own comment)."""
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.excitation_energies/excitation_gap: "
                "only itensor_version=\"python\" is supported -- see "
                "pyitensor.idmrg_excitations' own module docstring")
        if self._result is None:
            self.gs_energy()
        if self._excitation_env is None:
            from .pyitensor import idmrg_excitations
            self._excitation_env = idmrg_excitations.build_excitation_environment(
                self._result, self._h_intra.op, self._h_inter.op, self.n_uc,
                self.site_types)
        return self._excitation_env

    def excitation_energies(self, k, n=1):
        """The lowest `n` excitation energies (above the ground state) at
        momentum `k` (radians, per unit cell) of the tangent-space/
        quasiparticle excitation ansatz -- see
        pyitensor.idmrg_excitations' own module docstring for the algorithm,
        and its "KNOWN LIMITATION" section for an important scope
        restriction: only a product-state-like (bond dimension D=1)
        converged ground state is currently supported (raises
        NotImplementedError otherwise) -- a genuinely entangled ground state
        was found, during development, to give a dispersion that is
        anomalously flat compared to the expected answer, not yet
        root-caused."""
        env = self._get_excitation_environment()
        from .pyitensor import idmrg_excitations
        return idmrg_excitations.excitation_energies(env, k, n=n)

    def excitation_gap(self, ks=None):
        """The scalar "first excited state" gap: min_k
        excitation_energies(k, n=1) over a momentum scan (default
        `ks=numpy.linspace(-pi, pi, 41)`) -- mirrors the finite-chain
        `get_gap()` naming (docs/user_guide.md Sec.4), which returns
        `get_excited(n=2)[1]-get_excited(n=2)[0]`; here there is no
        discrete second state, only a momentum-resolved band, so the gap
        is defined as the band's own minimum. See excitation_energies'
        own docstring for the scope restriction (D=1 ground states only,
        for now)."""
        if ks is None:
            ks = np.linspace(-np.pi, np.pi, 41)
        return min(self.excitation_energies(k, n=1)[0] for k in ks)

    def local_excitation_gap(self, window=0, niter=200):
        """A cheap, cruder alternative to excitation_gap: re-diagonalizes
        the growing algorithm's own final 2-site effective Hamiltonian
        (already solved once for the ground state) for its second-lowest
        eigenvalue, and returns the difference -- the "local superblock
        gap". Unlike excitation_gap (the tangent-space/quasiparticle
        ansatz), this has no momentum label, does not require D=1, and
        reuses the ground state's own HL/HR environments unmodified rather
        than letting them relax for the excited sector -- see
        pyitensor.idmrg.local_excitation_gap's own docstring for the full
        rationale (including why this, not a soft penalty weight, is the
        exact analogue of finite-chain DMRG's Lagrange-multiplier excited-
        state trick here) and for an important accuracy caveat: it is not
        guaranteed to match the true minimum-momentum gap the way
        excitation_gap does.

        `window` (default 0, the original behavior above) grows the local
        diagonalization block by `window` extra *free* physical sites on
        each side of the original 2 (both the ground state and the deflated
        first excited state are then re-solved fresh within this larger
        block, rather than reusing the growing algorithm's own ground
        vector) -- see pyitensor.idmrg.local_excitation_gap_windowed's own
        docstring for the construction. This measurably tightens the gap:
        on a D=1 field-polarized XX chain (exact answer known) the error
        drops from 10% at window=0 to <1% by window=4-6; on a genuinely
        entangled (D>1) transverse-field Ising chain the window=3 estimate
        matches an 18-site open finite chain's own ED gap to <1%, converging
        at least as fast as growing the finite chain itself does. Costs
        d**(2*window) more in local Hilbert space dimension per extra site
        pair (d = the physical dimension), so window>2-3 gets expensive
        quickly for d>2 (e.g. S=1). Only `n_uc=1` is supported for
        `window>0` (raises NotImplementedError otherwise -- see
        pyitensor.idmrg.local_excitation_gap_windowed's own docstring for
        why n_uc=2 needs sublattice-position bookkeeping not implemented
        yet).

        Only supported for itensor_version="python" -- see vev's own
        comment."""
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.local_excitation_gap: only "
                "itensor_version=\"python\" is supported -- see "
                "pyitensor.idmrg.local_excitation_gap's own docstring")
        if self._result is None:
            self.gs_energy()
        from .pyitensor import idmrg
        if window == 0:
            return idmrg.local_excitation_gap(self._result, niter=niter)
        return idmrg.local_excitation_gap_windowed(
            self._result, self._h_intra.op, self._h_inter.op, self.site_types,
            self.n_uc, window=window, niter=niter)

    def kpm_finite(self, opname_i, p_i, opname_j, r, n_window,
                   window_chain_kwargs=None, **kwargs):
        """Dynamical correlator <opname_i(site p_i) opname_j(site p_i+r)>
        (omega) of the infinite chain, computed with the KPM method --
        reusing the existing finite-chain KPM implementation
        (kpmdmrg.get_dynamical_correlator, pyitensor/chain.py's
        Chain.kpm_dynamical_correlator) verbatim, not a new Chebyshev-
        recursion implementation. Named `kpm_finite` (rather than
        `get_dynamical_correlator`, the finite-chain method it wraps) to
        flag up front that this is the finite-window *approximation*
        described below, not an exact infinite-size calculation -- see
        "A genuinely infinite-chain KPM (not implemented here)" at the end
        of this docstring for what an exact version would need.

        == Method: finite window, open boundary conditions ==

        vev/correlator work directly with the converged, exactly
        translationally-invariant unit cell (IDMRGResult.U_list) via the
        transfer-matrix formalism -- but there is no analogous *dynamical*
        formalism here: H is extensive/unbounded in the thermodynamic
        limit, so a literal Chebyshev expansion of the full infinite H has
        no meaning (unlike apply_mpo/imps_sum's *bounded*-operator scope).
        Instead, this method builds an ordinary finite, open-boundary chain
        of `n_window` repeats of this chain's own unit cell
        (n_window*n_uc physical sites total, Hamiltonian = h_intra tiled
        onto every cell plus h_inter tiled onto every adjacent pair of
        cells -- see _window_hamiltonian), places opname_i/opname_j at
        sites p_i/p_i+r of the window's *central* unit cell (as far as
        possible from both open ends), and delegates directly to
        kpmdmrg.get_dynamical_correlator on that ordinary finite
        Many_Body_Chain (itensor_version="python") -- the exact same KPM
        code path an ordinary finite Spin_Chain/Fermionic_Chain uses
        today. (kpmdmrg.get_dynamical_correlator is called directly,
        bypassing dynamics.py's own is_hermitian() gate -- see
        set_hamiltonian's own comment for why that check is known to
        false-reject an ordinary cross-site Heisenberg-style term.)

        This is a genuine approximation to the (only ill-defined via KPM
        anyway) infinite-chain dynamical correlator, not an exact
        infinite-size method: results carry finite-size/open-boundary
        corrections that must be checked by convergence in `n_window`,
        exactly as a static vev/correlator caller would check maxm/etol
        convergence of the original iDMRG ground state. In particular: one
        Chebyshev moment corresponds to one application of the (nearest-
        neighbor) window Hamiltonian, so it can only move information by
        ~1 site per moment (a Lieb-Robinson-style bound) -- but KPM's own
        moment count `n` (see kpmdmrg.get_dynamical_correlator /
        Chain._scaled_hamiltonian) scales with the *window's own extensive
        bandwidth* divided by the requested `delta` (an ordinary finite
        chain's KPM already has this property -- nothing new here), so a
        genuinely fine `delta` can require a moment count comparable to
        (or larger than) `n_window` itself, at which point boundary
        reflections contaminate the result regardless of how large
        `n_window` is -- there is no way around this for a fixed `delta`
        short of an actual infinite-boundary-condition (IBC) window
        method, out of scope here (see below). Prefer a coarser `delta`,
        or check that the correlator has visibly converged with growing
        `n_window`, for quantitative work (especially near a gapless
        point, where a fine `delta` is most tempting).

        opname_i/opname_j: named single-site operators, at sublattice
        position p_i (0..n_uc-1) of the window's central cell and p_i+r of
        whichever cell that lands in (r>=0 physical sites, see
        two_point_correlator's own r convention).

        n_window: number of unit-cell repeats in the finite window -- no
        default (see the convergence caveat above); must at least be able
        to fit p_i and p_i+r (checked below), but should be chosen much
        larger than that in practice.

        window_chain_kwargs: optional dict of attribute overrides applied
        to the temporary finite Many_Body_Chain before set_hamiltonian
        (e.g. dict(maxm=60, nsweeps=20, kpmmaxm=80, kpm_scale=0.5)) -- the
        same attributes an ordinary finite chain exposes
        (manybodychain.py's Many_Body_Chain.__init__), independent of this
        iDMRG chain's own self.maxm/etc.

        Remaining **kwargs (delta, kernel, es, deconvolve, ...) are
        forwarded to kpmdmrg.get_dynamical_correlator unchanged. Returns
        (es, correlator), exactly like a finite chain's own
        get_dynamical_correlator.

        Unlike vev/correlator, this does not need self._result (no
        dependency on a previously converged IDMRGResult, or even on
        self.itensor_version -- it only needs the Hamiltonian specification
        from set_hamiltonian), so it works regardless of which backend
        gs_energy() itself used.

        == A genuinely infinite-chain KPM (not implemented here) ==

        A KPM calculation that is exact in the thermodynamic limit (not
        just "checked to have converged in n_window") is possible in
        principle, via infinite boundary conditions (IBC, Phien/McCulloch/
        Vidal, PRB 86, 245107 (2012)) -- but is a substantially larger
        feature than this method, and is not implemented here:

        1. Cap the window's two ends with the *converged environment*
           (the HL/HR accumulated Hamiltonian blocks idmrg_ground_state
           already builds during growth, or an equivalent construction
           from the fixed-point machinery apply_mpo/imps_overlap already
           use) instead of plain open boundaries. This removes a real
           error source this method's open-boundary window has that
           `n_window` alone cannot fix no matter how large: an open
           chain's own ground state carries boundary artifacts (e.g.
           Friedel-oscillation-like features) that contaminate even the
           *central* region, not just the two edges, because the window's
           own DMRG ground-state search has no way to know it should
           match the true infinite bulk rather than terminate at a
           physical edge.
        2. Independently, removing the `n_window`-vs-moment-count
           constraint entirely (rather than just picking `n_window`
           comfortably larger than the expected moment count, as this
           method requires) needs the *active* window to grow by roughly
           one site per Chebyshev moment as the recursion proceeds
           (reusing the converged unit-cell tensors to extend the state
           on demand, similar in spirit to a "growing window" real-time
           TEBD/TDVP simulation on an infinite chain), rather than fixing
           a static `n_window` up front. Done this way, the calculation's
           cost is bounded only by the usual bond-dimension/computational
           budget of any DMRG calculation, not by a fixed geometric
           window guessed in advance.

        Both pieces reuse machinery this module already has (the HL/HR
        environment construction from idmrg_ground_state, and the
        fixed-point/transfer-matrix code apply_mpo/imps_overlap build
        on), but assembling them into a working IBC-window KPM
        implementation is a new, nontrivial feature in its own right, not
        a small patch on top of kpm_finite -- flagged as a known,
        deliberate follow-up rather than attempted here."""
        if self._h_intra is None:
            raise RuntimeError(
                "Infinite_Many_Body_Chain.kpm_finite called "
                "before set_hamiltonian")
        if not (0 <= p_i < self.n_uc):
            raise ValueError("kpm_finite: p_i must be in 0..{} "
                              "(n_uc-1), got {!r}".format(self.n_uc - 1, p_i))
        if r < 0:
            raise ValueError("kpm_finite: r must be >= 0")
        if n_window < 1:
            raise ValueError("kpm_finite: n_window must be >= 1")
        n_sites = n_window * self.n_uc
        center = n_window // 2
        s_i = center * self.n_uc + p_i
        s_j = s_i + r
        if s_j >= n_sites:
            raise ValueError(
                "kpm_finite: n_window={} is too small to fit "
                "opname_j at window site {} (p_i={}, r={}, central cell {}) "
                "-- increase n_window".format(n_window, s_j, p_i, r, center))
        from .manybodychain import Many_Body_Chain
        from . import kpmdmrg
        window_sites = self.site_types * n_window
        wc = Many_Body_Chain(window_sites, itensor_version="python")
        for k, v in (window_chain_kwargs or {}).items():
            setattr(wc, k, v)
        h_window = _window_hamiltonian(self._h_intra, self._h_inter, self.n_uc, n_window)
        wc.set_hamiltonian(h_window)
        op_i = wc.get_operator(opname_i, s_i)
        op_j = wc.get_operator(opname_j, s_j)
        return kpmdmrg.get_dynamical_correlator(wc, name=(op_i, op_j), **kwargs)

    def td_dynamical_correlator(self, opname_i, p_i, opname_j, n_window,
                                 dt=0.1, nt=200, x_values=None,
                                 maxdim=60, cutoff=1e-10, niter=50,
                                 connected=True, **kwargs):
        """Real-time dynamical correlator `S(k,omega)` of the infinite
        chain, via infinite-boundary-condition (IBC) window time evolution
        (Milsted/Vanderstraeten et al., "Infinite boundary conditions for
        response functions and limit cycles in iDMRG", arXiv:1804.09163,
        Sec. V.1) -- computed for *every* k and every distance x from a
        *single* window TDVP evolution (`pyitensor.idmrg_window`'s own
        `dynamical_correlator_komega`), unlike `kpm_finite` (one run per
        `(opname_i, opname_j, r)` triple).

        This is the actual infinite-boundary-condition method
        `kpm_finite`'s own docstring flags as "not implemented here" (its
        own "A genuinely infinite-chain KPM" section, item 1: cap the
        window's two ends with the converged environment instead of plain
        open boundaries) -- built on real-time TDVP rather than
        KPM/Chebyshev, so item 2 of that same section ("growing window",
        needed there to remove a moment-count-vs-`n_window` constraint) is
        naturally moot here too: the window is IBC-capped with the *true*
        converged environment from the start, so there is no open-boundary
        artifact (e.g. Friedel-oscillation-like contamination of the
        window's own central region) to grow away from in the first
        place, unlike `kpm_finite`'s own open-boundary window.

        `itensor_version="python"` (needs `pyitensor.idmrg`'s own
        `IDMRGResult.env_HL`/`env_HR` snapshot) and `itensor_version=3`
        (native `mpscpp3.Chain::td_dynamical_correlator_window`, reusing
        ITensorTDVP's own boundary-tensor `tdvp()` overload directly rather
        than pyitensor's hand-rolled window-aware sweep -- see that
        method's own comment) are both supported; any other backend raises
        NotImplementedError. Calls `gs_energy()` first if not already run
        (`self._result`/`self._session3` unset, per backend). The v3 path's
        own `x_values` may not extend beyond the window's own explicit
        range (`center+x` must stay within `[1, n_window*n_uc]`) -- unlike
        the "python" backend, it does not pad beyond the window with extra
        unevolved unit-cell copies, so increase `n_window` instead if a
        wider `x_values` is needed.

        `opname_j` is applied at sublattice position `p_i` (0..n_uc-1) and
        evolved forward in time; `opname_i` is inserted at the shifted
        position `p_i`'s own site `+x` (bra side, *not* time-evolved) --
        see `pyitensor.idmrg_window.dynamical_correlator_td`'s own
        docstring for the precise (Schrödinger-picture, background-
        subtracted) convention this follows, and why no `e^{iE0t}`
        Heisenberg-picture conversion is applied (matching this codebase's
        own established "TD" submode convention).

        `n_window`: number of unit-cell repeats in the window -- same
        convergence caveat as `kpm_finite`'s own `n_window` (check by
        increasing it and confirming results stop changing), except here
        the two ends are physically correct (IBC-capped) rather than
        open-boundary artifacts, so a much smaller margin should suffice
        in practice. `dt`/`nt`: TDVP time step and step count (total
        simulated time `nt*dt`); `maxdim`/`cutoff`/`niter`: TDVP
        truncation/Krylov-dimension controls (see
        `pyitensor.idmrg_window.window_tdvp_step`). `x_values`: distances
        (physical sites, relative to `opname_j`'s own site) to reconstruct
        `S(x,t)` for before Fourier transforming -- defaults to a
        conservative margin based on `n_window` (see
        `dynamical_correlator_td`'s own docstring); must stay within the
        causal cone for `nt*dt`, same spirit as `n_window` itself.
        `connected=True` (default) subtracts the disconnected background
        `<opname_i><opname_j>` -- confirmed directly to matter (see
        `dynamical_correlator_td`'s own docstring): without it, the
        spatial Fourier transform is dominated by a spurious k=0
        contribution with no discernible dispersion structure.

        Remaining `**kwargs` (`ks`, `es`, `delta`, `window`, `factor`) are
        forwarded to `pyitensor.idmrg_window.dynamical_correlator_komega`
        unchanged -- `delta` means the same Lorentzian-broadening thing
        here as in every other dynamical-correlator submode in this
        codebase (it reuses `timedependent.py`'s own
        `_fourier_transform_correlator` directly). Returns `(ks, es,
        Skw)`, `Skw` shaped `(len(ks), len(es))`."""
        if not (0 <= p_i < self.n_uc):
            raise ValueError("td_dynamical_correlator: p_i must be in 0..{} "
                              "(n_uc-1), got {!r}".format(self.n_uc - 1, p_i))
        if self.itensor_version == "python":
            if self._result is None:
                self.gs_energy()
            from .pyitensor import idmrg_window
            return idmrg_window.dynamical_correlator_komega(
                self._result, n_window, opname_i, opname_j, dt, nt,
                cutoff=cutoff, maxdim=maxdim, niter=niter, x_values=x_values,
                connected=connected, p_i=p_i, **kwargs)
        if self.itensor_version != 3:
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.td_dynamical_correlator: only "
                "itensor_version=\"python\" or itensor_version=3 are "
                "supported")
        if self._session3 is None:
            self.gs_energy()
        # Same default margin idmrg_window.py's own dynamical_correlator_td
        # uses -- computed here since the v3 binding needs a concrete list
        # (no Python-side IDMRGResult to read n_uc off of on this path).
        if x_values is None:
            margin = max(1, n_window * self.n_uc // 4)
            x_values = range(-margin, margin + 1)
        xs_in = sorted(x_values)
        ts, xs, S = self._session3.td_dynamical_correlator_window(
            n_window, opname_i, opname_j, dt, nt, xs_in, maxdim, cutoff,
            niter, connected, p_i)
        ks = kwargs.pop("ks", None)
        es = kwargs.pop("es", None)
        delta = kwargs.pop("delta", 5e-2)
        window = kwargs.pop("window", [-1, 10])
        factor = kwargs.pop("factor", 1)
        if kwargs:
            raise TypeError("td_dynamical_correlator: unexpected kwargs {!r} "
                             "for itensor_version=3".format(list(kwargs)))
        from .timedependent import _fourier_transform_correlator
        if ks is None:
            ks = np.linspace(-np.pi, np.pi, 200)
        ks = np.asarray(ks)
        Skw = None
        es_out = es
        for ik, k in enumerate(ks):
            phase = np.exp(-1j * k * xs)
            Skt = S @ phase
            es_k, gk = _fourier_transform_correlator(
                ts, Skt, dt, es=es_out, window=window, delta=delta, factor=factor)
            if Skw is None:
                es_out = es_k
                Skw = np.zeros((len(ks), len(es_k)), dtype=complex)
            Skw[ik] = gk
        return ks, es_out, Skw


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
