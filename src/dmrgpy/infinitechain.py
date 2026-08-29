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
something else. `self.gs_method` (`"vumps"` by default, on EITHER
backend -- see `Infinite_Many_Body_Chain`'s own class docstring) picks
the ground-state algorithm; `vev`/`correlator` work for BOTH values on
BOTH backends: `"vumps"` via `pyitensor.vumps.onsite_expectation`/
`two_point_correlator` and `Chain::vumps_onsite_expectation`/
`vumps_two_point_correlator` (the AC/AR-based formula), `"idmrg"` via
`pyitensor.idmrg`'s and `Chain::idmrg_onsite_expectation`/
`idmrg_two_point_correlator`'s own dominant-fixed-point-based versions
over the gauge-consistent unit cell. `local_excitation_gap` needs
`gs_method="idmrg"` and works on both backends at `window=0`
(`window>0` stays `"python"`-only, see that method's own docstring);
`excitation_energies`/`excitation_gap` need `gs_method="vumps"` and work
on both backends. See `vev`'s own docstring for the full support
matrix.

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

import warnings

import numpy as np

from . import multioperator
from .multioperator import MultiOperator
from .pyitensor.sites.base import is_fermionic


def _warn_if_growth_missed_local_ground_state(e0_fresh, e0_stored):
    """`local_excitation_gap` re-solves the stored 2-site effective
    Hamiltonian for its own ground state rather than reusing the growing
    algorithm's reported one (see
    `pyitensor.idmrg._lowest_two_eigenvalues`'s docstring for why: that
    loop's convergence test bounds the distance to *some* eigenvalue, not
    to the lowest one, so its number and a strong random-started solve here
    need not be eigenvalues of the same eigenpair). The two agree whenever
    the growth loop did find that Hamiltonian's ground state, which is the
    ordinary case.

    When they do not, the returned gap is still a genuine spectral gap of
    the stored operator, but it is no longer a gap measured above the state
    whose observables `vev`/`correlator` report: the growth loop warm-starts
    every micro-step from the previous one and so stays on the branch it
    started on, while a random-started solve of the same operator is free to
    leave it -- and with `ConserveQNs=false` (both backends) nothing confines
    that search to the converged state's own particle-number sector, so what
    it finds can be a charge-imbalanced state against the frozen environment.
    Worth a warning rather than a silent number."""
    if abs(e0_fresh - e0_stored) > 1e-6 * (1.0 + abs(e0_stored)):
        warnings.warn(
            "local_excitation_gap: the growing algorithm's own final local "
            "ground eigenvalue ({:.10g}) is above this operator's actual "
            "ground eigenvalue ({:.10g}), by {:.4g}. The returned gap is the "
            "true spectral gap of that stored 2-site effective Hamiltonian, "
            "but it is NOT a gap above the converged state vev()/correlator() "
            "describe -- with ConserveQNs=false the local spectrum also holds "
            "charge-imbalanced states against the frozen environment, which "
            "the warm-started growth loop never visits. Treat this gap as "
            "unreliable for this model.".format(
                e0_stored, e0_fresh, e0_stored - e0_fresh),
            RuntimeWarning, stacklevel=3)


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
        # Odd total fermion parity: the term's Jordan-Wigner string opens
        # and never closes, so it would have to run to infinity in both
        # tiling directions -- unrepresentable by either backend's
        # per-term string threading (pyitensor/idmrg.py's
        # _term_site_matrices, mpscpp3/chain_session.h's
        # idmrg_classify_terms), and not a physical (parity-conserving)
        # Hamiltonian term anyway. Both backends reject it themselves, but
        # only after set_hamiltonian has returned: the v3 one does it with
        # an ITensor Error(), which aborts the whole process rather than
        # raising into Python. Checking it here, on the raw terms every
        # backend shares, turns that into an ordinary ValueError at the
        # point where the offending term is still in front of the user.
        if sum(1 for name, _site in term[1:] if is_fermionic(name)) % 2:
            raise ValueError(
                "Infinite_Many_Body_Chain.set_hamiltonian: a term has odd "
                "total fermion parity ({}) -- its Jordan-Wigner string "
                "cannot be closed within the unit cell".format(term))
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
    backend (`mpscpp3/chain_session.h`'s `Chain::idmrg_ground_state`/
    `Chain::vumps_ground_state`) computes the energy density under either
    `gs_method`, and supports `vev`/`correlator` under either too
    (`Chain::vumps_onsite_expectation`/`vumps_two_point_correlator`, and
    `Chain::idmrg_onsite_expectation`/`idmrg_two_point_correlator`).
    `itensor_version=2` and `"julia_live"`/`"julia"` have no iDMRG port at
    all and raise `NotImplementedError`.

    `self.gs_method`: `"vumps"` (default -- `pyitensor/vumps.py`'s/
    `Chain::vumps_ground_state`'s direct fixed-bond-dimension solver, see
    `vumps.py`'s own module docstring for the algorithm) or `"idmrg"`
    (`pyitensor/idmrg.py`'s/`Chain::idmrg_ground_state`'s growing/
    infinite-size algorithm). `vev`/`correlator` work under EITHER
    `gs_method` on `itensor_version="python"` -- `"idmrg"` dispatches to
    `pyitensor.idmrg.onsite_expectation`/`two_point_correlator` (via
    `IDMRGResult.cell_list`, the gauge-consistent unit cell, and a
    dominant-fixed-point eigenproblem), `"vumps"` to
    `pyitensor.vumps`'s own same-named functions (via `VUMPSResult`'s
    mixed-gauge `{AC, AR}` directly, no eigenproblem needed -- see those
    functions' own docstrings) -- and under EITHER `gs_method` on
    `itensor_version=3` as well, via the C++ ports of those same two
    families. `local_excitation_gap` AND `td_dynamical_correlator` both
    require `gs_method="idmrg"` specifically, on EITHER backend (the
    former re-diagonalizes the growing algorithm's own final 2-site
    effective Hamiltonian; the latter needs the growing algorithm's own
    converged environment snapshot -- neither has a VUMPS equivalent);
    `local_excitation_gap`'s `window>0` variant is additionally
    `itensor_version="python"`-only, see its own docstring;
    conversely `excitation_energies`/`excitation_gap` (the tangent-space/
    quasiparticle excitation ansatz) require `gs_method="vumps"` -- they
    need `VUMPSResult`'s own mixed-gauge {AL,AR,C,GL,GR}, which the
    growing algorithm's `IDMRGResult` has no equivalent of (see
    `pyitensor.idmrg_excitations`' own module docstring). Since `"vumps"`
    is the default, code that specifically needs the growing algorithm
    (`local_excitation_gap`, `td_dynamical_correlator`, or just a
    preference for `idmrg_ground_state`'s more battle-tested D>1
    convergence -- see `gs_energy`'s own reliability comment) MUST set
    `self.gs_method = "idmrg"` explicitly; it is no longer implied by
    doing nothing."""

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
        # n_uc > 2 is supported through gs_method="vumps" (the default),
        # which uses the sequential multi-site algorithm -- see
        # pyitensor/vumps_ms.py. It is NOT supported by gs_method="idmrg",
        # whose growth loop pairs sublattice p_L=mstep with
        # p_R=n_uc-1-mstep and so only produces genuinely adjacent active
        # sites for n_uc in {1,2}; for n_uc>=3 it pairs non-adjacent
        # positions whose automaton channels carry different, bond-specific
        # pending content, which would silently mix unrelated bond content
        # rather than fail to contract. That path rejects n_uc>2 itself
        # (gs_energy below, and pyitensor/idmrg.py's own guard), so the
        # restriction now lives with the algorithm that actually has it
        # instead of blocking construction outright.

        self.gs_method = "vumps"  # "vumps" (the default -- pyitensor.vumps's/
                                  # Chain::vumps_ground_state's direct
                                  # fixed-bond-dimension solver, supported
                                  # on both itensor_version="python" and 3)
                                  # or "idmrg" (pyitensor.idmrg's/
                                  # Chain::idmrg_ground_state's growing/
                                  # infinite-size algorithm) -- see
                                  # gs_energy's own comment. local_excitation_
                                  # gap/td_dynamical_correlator need
                                  # "idmrg" specifically; excitation_energies/
                                  # excitation_gap need "vumps" specifically
                                  # -- see this class's own docstring for the
                                  # full support matrix. Mirrors
                                  # manybodychain.py's Many_Body_Chain.
                                  # tevol_method pattern (a plain attribute
                                  # picking between two algorithms for the
                                  # same public method).
        self.maxm = 30          # wavefunction bond dimension cap
        self.cutoff = 1e-12     # SVD truncation
        self.maxiter = 200      # iDMRG macro-iterations (growth steps);
                                 # gs_method="vumps": VUMPS outer-iteration
                                 # cap per bond dimension in its own D-ramp
                                 # (pyitensor.vumps.vumps_ground_state's own
                                 # `maxiter`)
        self.etol = 1e-10       # energy-density convergence tolerance;
                                 # gs_method="vumps": VUMPS's own
                                 # gauge-mismatch convergence tolerance
                                 # (`tol`).
                                 # gs_method="idmrg": meaningful down to
                                 # machine precision since idmrg.py started
                                 # subtracting a per-site energy baseline
                                 # from both growing environments
                                 # (`_subtract_energy_baseline`, the
                                 # equivalent of its ITensor reference's
                                 # `HL += -energy*IL`). Without it the
                                 # superblock eigenvalue grew extensively
                                 # while _lanczos_ground_state stopped on a
                                 # *relative* criterion, so the absolute
                                 # noise in the finite-difference density
                                 # grew linearly with the iteration count:
                                 # measured on a gapped n_uc=2 chain,
                                 # |E| = 603 after 400 macro-iterations with
                                 # the density jittering over ~9e-11 long
                                 # after it had physically converged -- right
                                 # at this etol. It is now |E| = 3.8 and
                                 # ~4e-15.
        self.vumps_nrestarts = 4  # gs_method="vumps" only: independent
                                   # random-restart attempts per bond
                                   # dimension in VUMPS's own D-ramp -- see
                                   # pyitensor.vumps.vumps_ground_state's
                                   # own docstring for why this matters
                                   # (single-attempt VUMPS from a random
                                   # start is not reliable for D>1)
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
        self.noise = 1e-4       # gs_method="idmrg" only: White density-matrix
                                 # noise, applied ONLY while the growing
                                 # state is still a product state (plus a
                                 # short tail), so an already-entangled model
                                 # never sees it and is unaffected. Exists
                                 # because a product state is an EXACT fixed
                                 # point of the growth loop that no amount of
                                 # solver refinement escapes -- see
                                 # pyitensor.idmrg._noise_perturbed_split and
                                 # docs/known_issue_idmrg_product_state_
                                 # collapse.md. 0 disables it entirely (and
                                 # restores the pre-fix numerics exactly).
        self.noise_iters = 40   # cap on how many macro-iterations may carry
                                 # noise, so a model whose ground state
                                 # genuinely IS a product state (a polarized
                                 # chain) stops re-arming it and still ends
                                 # on a clean, unperturbed tail.
        self.verbose = False

        self.hamiltonian = None  # user-facing MultiOperator (L/C/R indices)
        self._h_intra = None
        self._h_inter = None
        self._result = None      # pyitensor.idmrg.IDMRGResult once converged
                                  # (itensor_version="python", gs_method=
                                  # "idmrg" only -- itensor_version=3 keeps
                                  # the equivalent state inside its own C++
                                  # Chain instead, see self._session3 below,
                                  # and supports vev/correlator from it just
                                  # the same)
        self._vumps_result = None  # pyitensor.vumps.VUMPSResult once
                                    # converged (itensor_version="python",
                                    # gs_method="vumps" only -- v3's own
                                    # gs_method="vumps" result lives inside
                                    # self._session3 instead, see below)
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
                                  # Chain instance gs_energy() last ran a
                                  # ground-state method on, kept alive so
                                  # every observable can reuse its private
                                  # snapshots -- after idmrg_ground_state:
                                  # the converged environments
                                  # (td_dynamical_correlator), the
                                  # gauge-consistent unit cell (vev/
                                  # correlator) and the final local
                                  # superblock (local_excitation_gap);
                                  # after vumps_ground_state: the
                                  # mixed-gauge state (vev/correlator/
                                  # excitation_energies). None of these is
                                  # exposed to Python directly -- see
                                  # chain_session.h's own comments -- so a
                                  # *fresh* Chain() has none of them, and
                                  # re-running gs_energy() on a new
                                  # instance would silently lose them.
        self._session3_has_vumps = False  # itensor_version=3 only: whether
                                           # self._session3 last ran
                                           # vumps_ground_state (True) or
                                           # idmrg_ground_state (False) --
                                           # excitation_energies' own v3
                                           # branch needs this (not just
                                           # "is self._session3 None") to
                                           # correctly rerun gs_energy() if
                                           # gs_method was switched to
                                           # "vumps" after an earlier
                                           # gs_method="idmrg" run already
                                           # populated self._session3 --
                                           # mirrors self._vumps_result's
                                           # own None-sentinel role for the
                                           # "python" backend above.

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
        self._vumps_result = None
        self.e0 = None
        self.converged = None
        self._excitation_env = None
        self._session3 = None
        self._session3_has_vumps = False

    def gs_energy(self):
        """Run the ground-state solver (self.maxiter iterations/macro-
        iterations at most) and return the ground-state energy *per site*
        -- the physically meaningful observable for an infinite system (a
        total energy would be unboundedly large).

        itensor_version="python" dispatches on `self.gs_method`:
        "vumps" (the default) runs pyitensor/vumps.py's direct fixed-bond-
        dimension solver (see its own module docstring for the algorithm
        -- in particular, its `.converged` diagnostic and general
        robustness at larger D are both less battle-tested than
        idmrg_ground_state's, a real tradeoff for defaulting to it, see
        that module's own "Convergence robustness" section), keeping the
        resulting VUMPSResult in self._vumps_result; "idmrg" instead runs
        pyitensor/idmrg.py's own growing algorithm and keeps the resulting
        IDMRGResult in self._result for vev/correlator to reuse (see those
        methods) -- required for local_excitation_gap/
        td_dynamical_correlator specifically, which have no VUMPS
        equivalent (see this class's own docstring). vev/correlator
        support BOTH gs_method values here, computed directly from AC/AR
        for "vumps" (no per-sublattice U_list needed, see vev's own
        docstring). itensor_version=3 now ALSO dispatches on
        `self.gs_method`, mirroring the "python" backend: "vumps" (default)
        calls the compiled mpscpp3 backend's Chain::vumps_ground_state
        (C++ port of pyitensor/vumps.py -- see mpscpp3/chain_session.h's
        own doc comment at that method for the algorithm/scope, including
        the one simplification it takes relative to pyitensor's own
        driver); "idmrg" calls Chain::idmrg_ground_state instead. Both
        support vev/correlator on this backend -- "vumps" via
        Chain::vumps_onsite_expectation/vumps_two_point_correlator (a C++
        port of pyitensor/vumps.py's own AC/AR-based formula), "idmrg" via
        Chain::idmrg_onsite_expectation/idmrg_two_point_correlator (a C++
        port of pyitensor/idmrg.py's own transfer-matrix machinery over
        the gauge-consistent unit cell) -- see vev's own docstring for the
        full support matrix. self._result/self._vumps_result are both left
        None on the itensor_version=3 path either way (they are
        "python"-backend-only caches): everything the v3 backend's own
        observables need lives inside the C++ Chain, as private snapshots
        set by whichever ground-state method ran. self._session3 is
        therefore kept alive either way -- td_dynamical_correlator and
        local_excitation_gap need it after gs_method="idmrg",
        excitation_energies/excitation_gap after gs_method="vumps", and
        vev/correlator after either (see _get_excitation_environment's own
        docstring)."""
        if self._h_intra is None:
            raise RuntimeError(
                "Infinite_Many_Body_Chain.gs_energy called before set_hamiltonian")
        self._excitation_env = None  # about to rebuild self._result; see __init__'s own comment
        if self.itensor_version == "python" and self.gs_method == "idmrg":
            from .pyitensor import idmrg
            self._result = idmrg.idmrg_ground_state(
                self.site_types, self._h_intra.op, self._h_inter.op, self.n_uc,
                maxm=self.maxm, cutoff=self.cutoff, maxiter=self.maxiter,
                etol=self.etol, niter=self.niter, verbose=self.verbose,
                noise=self.noise, noise_iters=self.noise_iters)
            self._vumps_result = None
            self.e0 = self._result.e0
            self.converged = self._result.converged
        elif self.itensor_version == "python" and self.gs_method == "vumps":
            from .pyitensor import vumps
            self._vumps_result = vumps.vumps_ground_state(
                self.site_types, self._h_intra.op, self._h_inter.op, self.n_uc,
                D=self.maxm, tol=self.etol, maxiter=self.maxiter,
                niter_lanczos=self.niter, nrestarts=self.vumps_nrestarts,
                verbose=self.verbose)
            self._result = None
            self.e0 = self._vumps_result.e0
            self.converged = self._vumps_result.converged
        elif self.itensor_version == "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.gs_energy: gs_method={!r} is not "
                "implemented -- only \"idmrg\" and \"vumps\" are supported "
                "for itensor_version=\"python\"".format(self.gs_method))
        elif self.itensor_version == 3 and self.gs_method == "idmrg":
            chain = self._make_cpp_chain()
            # jordan_wigner_transform=False: the C++ backend threads the
            # Jordan-Wigner string itself (Chain::idmrg_classify_terms, a
            # port of pyitensor/idmrg.py's _term_site_matrices), locally
            # between each term's own two endpoints, exactly like the
            # itensor_version="python" branches above -- see
            # MultiOperator.to_terms' own docstring for why the finite-chain
            # transform is wrong here.
            terms_intra = self._h_intra.to_terms(jordan_wigner_transform=False)
            terms_inter = self._h_inter.to_terms(jordan_wigner_transform=False)
            density, converged, _niter_done = chain.idmrg_ground_state(
                terms_intra, terms_inter, self.maxm, self.cutoff,
                self.maxiter, self.etol, self.niter, self.restarts,
                self.noise, self.noise_iters)
            self._session3 = chain  # kept alive for td_dynamical_correlator -- see __init__'s own comment
            self._session3_has_vumps = False
            self._result = None
            self._vumps_result = None
            self.e0 = density
            self.converged = converged
        elif self.itensor_version == 3 and self.gs_method == "vumps":
            chain = self._make_cpp_chain()
            # jordan_wigner_transform=False -- see the "idmrg" branch above
            # (vumps_ground_state classifies its terms with the very same
            # Chain::idmrg_classify_terms).
            terms_intra = self._h_intra.to_terms(jordan_wigner_transform=False)
            terms_inter = self._h_inter.to_terms(jordan_wigner_transform=False)
            # self.niter feeds the H_AC/H_C Krylov solves here exactly as
            # it does pyitensor.vumps.vumps_ground_state's own
            # niter_lanczos; it only bites once those solves are big enough
            # to go matrix-free (Chain::vx_dense_eig_max_).
            e0, converged, _niter_done, _gauge_mismatch = chain.vumps_ground_state(
                terms_intra, terms_inter, self.maxm, self.etol, self.maxiter,
                self.vumps_nrestarts, max(self.niter, 2))
            self._session3 = chain  # kept alive for excitation_energies/excitation_gap
            self._session3_has_vumps = True
            self._result = None
            self._vumps_result = None
            self.e0 = e0
            self.converged = converged
        else:  # itensor_version == 3, unknown gs_method
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.gs_energy: gs_method={!r} is not "
                "implemented -- only \"idmrg\" and \"vumps\" are supported "
                "for itensor_version=3".format(self.gs_method))
        return self.e0

    def _make_cpp_chain(self):
        """A fresh mpscpp3 Chain for this unit cell -- shared setup
        (backend-availability check + verbosity) for gs_energy's own
        itensor_version=3 branches (gs_method="idmrg"/"vumps")."""
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
        return chain

    def vev(self, opname, p, group="C"):
        """<opname> at site p (0..n_uc-1) of cell-group `group` -- by
        translational invariance this is identical for every group, `group`
        is accepted purely so callers can write whichever reads most
        naturally (SxL/SxC/SxR all describe the same infinite chain).

        All four (backend, gs_method) combinations are supported:
        gs_method="idmrg" dispatches to pyitensor.idmrg.onsite_expectation
        for itensor_version="python" and to
        Chain::idmrg_onsite_expectation for itensor_version=3;
        gs_method="vumps" to pyitensor.vumps.onsite_expectation and
        Chain::vumps_onsite_expectation respectively (via VUMPSResult's
        mixed-gauge AC directly -- see those functions' own docstrings).

        gs_method="idmrg"'s static observables tile the gauge-consistent
        unit cell extracted from a single micro-step's theta
        (pyitensor/idmrg.py's `_theta_cell` / IDMRGResult.cell_list, and
        Chain::idmrg_theta_cell on the C++ side), NOT the raw per-
        micro-step per-sublattice factors. Tiling those -- which is what
        this did before that extraction existed -- put <Sz> at -0.13
        against an exact 0 and <Sz(0)Sz(1)> at +0.028 against an exact
        -0.101 (the wrong sign) on the XX chain; it now reproduces both to
        the bond dimension's own truncation error."""
        if group not in ("L", "C", "R"):
            raise ValueError("vev: group must be 'L', 'C' or 'R', got {!r}".format(group))
        if not (0 <= p < self.n_uc):
            raise ValueError("vev: p must be in 0..{} (n_uc-1), got {!r}".format(
                self.n_uc - 1, p))
        if self.itensor_version == "python" and self.gs_method == "idmrg":
            if self._result is None:
                self.gs_energy()
            from .pyitensor import idmrg
            return idmrg.onsite_expectation(self._result, opname, p)
        if self.itensor_version == "python" and self.gs_method == "vumps":
            if self._vumps_result is None:
                self.gs_energy()
            from .pyitensor import vumps
            return vumps.onsite_expectation(self._vumps_result, opname, p)
        if self.itensor_version == 3 and self.gs_method == "vumps":
            if self._session3 is None or not self._session3_has_vumps:
                self.gs_energy()
            # n_uc>2 keeps a per-site (multi-site) snapshot rather than the
            # grouped one -- different representations, different reader.
            if self.n_uc > 2:
                return self._session3.vms_onsite_expectation(opname, p)
            return self._session3.vumps_onsite_expectation(opname, p)
        if self.itensor_version == 3 and self.gs_method == "idmrg":
            if self._session3 is None or self._session3_has_vumps:
                self.gs_energy()
            return self._session3.idmrg_onsite_expectation(opname, p)
        raise NotImplementedError(
            "Infinite_Many_Body_Chain.vev: only itensor_version=\"python\" "
            "and itensor_version=3, each with gs_method \"idmrg\" or "
            "\"vumps\", support static correlators -- got "
            "itensor_version={!r}, gs_method={!r}".format(
                self.itensor_version, self.gs_method))

    def correlator(self, opname_i, p_i, opname_j, r):
        """<opname_i(site p_i) opname_j(site p_i + r)>, r measured in
        physical sites (r>=0) -- see pyitensor.idmrg.two_point_correlator/
        pyitensor.vumps.two_point_correlator for the r=0 (same-site)
        convention.

        Same backend/gs_method support matrix as vev's own docstring
        (all four combinations of itensor_version in ("python", 3) with
        gs_method in ("idmrg", "vumps")), including its note on which
        tensors gs_method="idmrg" actually tiles.

        Fermionic operators: the Jordan-Wigner string is threaded across
        every site strictly between the two endpoints on every backend, so
        <Cdag_i C_j> here is the physical fermionic correlator, not a
        stringless product of two bare matrices."""
        if not (0 <= p_i < self.n_uc):
            raise ValueError("correlator: p_i must be in 0..{} (n_uc-1), got {!r}".format(
                self.n_uc - 1, p_i))
        if r < 0:
            # Checked here, not only per backend: the pyitensor paths raise
            # ValueError for this themselves, while the C++ ones raise an
            # ITError that pybind11 surfaces as a bare RuntimeError -- so
            # without this the same bad argument produced a different
            # exception type depending on itensor_version.
            raise ValueError(
                "Infinite_Many_Body_Chain.correlator: r must be >= 0 (use r "
                "and swap the operator order for the mirror separation), "
                "got {!r}".format(r))
        if (is_fermionic(opname_i) != is_fermionic(opname_j)):
            # Odd total fermion parity: the Jordan-Wigner string opened by
            # this pair is never closed, so it would have to run to infinity.
            # Such a correlator vanishes identically in any parity-conserving
            # state anyway; the backends reject it themselves (both C++
            # paths via Chain::idmrg_correlator_endpoints' own ITError), so
            # catch it here to give a clear ValueError on every backend.
            raise ValueError(
                "Infinite_Many_Body_Chain.correlator: the operator pair "
                "({!r}, {!r}) has odd total fermion parity -- its "
                "Jordan-Wigner string cannot be closed".format(
                    opname_i, opname_j))
        if self.itensor_version == "python" and self.gs_method == "idmrg":
            if self._result is None:
                self.gs_energy()
            from .pyitensor import idmrg
            return idmrg.two_point_correlator(self._result, opname_i, p_i, opname_j, r)
        if self.itensor_version == "python" and self.gs_method == "vumps":
            if self._vumps_result is None:
                self.gs_energy()
            from .pyitensor import vumps
            return vumps.two_point_correlator(self._vumps_result, opname_i, p_i, opname_j, r)
        if self.itensor_version == 3 and self.gs_method == "vumps":
            if self._session3 is None or not self._session3_has_vumps:
                self.gs_energy()
            if self.n_uc > 2:
                return self._session3.vms_two_point_correlator(
                    opname_i, p_i, opname_j, r)
            return self._session3.vumps_two_point_correlator(opname_i, p_i, opname_j, r)
        if self.itensor_version == 3 and self.gs_method == "idmrg":
            if self._session3 is None or self._session3_has_vumps:
                self.gs_energy()
            return self._session3.idmrg_two_point_correlator(
                opname_i, p_i, opname_j, r)
        raise NotImplementedError(
            "Infinite_Many_Body_Chain.correlator: only itensor_version="
            "\"python\" and itensor_version=3, each with gs_method "
            "\"idmrg\" or \"vumps\", support static correlators -- got "
            "itensor_version={!r}, gs_method={!r}".format(
                self.itensor_version, self.gs_method))

    def _get_excitation_environment(self):
        """Lazily build (or reuse) the pyitensor.idmrg_excitations.
        ExcitationEnvironment for the converged ground state -- expensive
        (a null-space computation plus several dense linear solves per
        momentum), so it is cached on self._excitation_env across repeated
        excitation_energies/excitation_gap calls (e.g. a k-scan), and
        invalidated whenever set_hamiltonian/gs_energy reruns (see
        __init__'s own comment).

        Requires itensor_version="python" with gs_method="vumps" (the
        default -- but still NOT satisfied automatically if a caller
        switched to gs_method="idmrg", e.g. for local_excitation_gap):
        the tangent-space excitation ansatz needs the
        mixed-gauge {AL,AR,C,GL,GR} representation only
        vumps.vumps_ground_state produces (see pyitensor.idmrg_excitations'
        own module docstring) -- the growing algorithm's IDMRGResult has no
        such representation. Set `self.gs_method = "vumps"` before calling
        gs_energy()/excitation_energies() (or let this method call
        gs_energy() itself, which will then also use gs_method="vumps").
        itensor_version=3 does NOT go through this method at all --
        excitation_energies/excitation_gap call the compiled mpscpp3
        backend's own Chain::vumps_excitation_energies directly instead
        (that backend keeps its own excitation-environment cache
        internally, mirroring this method's own caching one level down,
        see chain_session.h's own comment at
        vumps_build_excitation_environment)."""
        if self.n_uc > 2:
            raise NotImplementedError(
                "excitation_energies/excitation_gap: n_uc>2 (got {}) is not "
                "supported. The ground state itself is -- gs_method=\"vumps\" "
                "uses the sequential multi-site algorithm there "
                "(pyitensor/vumps_ms.py) -- but the tangent-space excitation "
                "ansatz built on top of it (pyitensor/idmrg_excitations.py) "
                "still expects the GROUPED single-supersite mixed gauge, so "
                "it would read a list of per-site tensors as one tensor. "
                "Porting the ansatz to the multi-site form is follow-up work; "
                "until then use local_excitation_gap (gs_method=\"idmrg\", "
                "n_uc<=2) or a smaller cell.".format(self.n_uc))
        if self.itensor_version != "python" or self.gs_method != "vumps":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.excitation_energies/excitation_gap: "
                "only itensor_version=\"python\" with gs_method=\"vumps\" is "
                "supported here -- set self.gs_method = \"vumps\" before "
                "calling gs_energy()/excitation_energies() -- see "
                "pyitensor.idmrg_excitations' own module docstring "
                "(itensor_version=3 is handled separately, directly inside "
                "excitation_energies/excitation_gap)")
        if self._vumps_result is None:
            self.gs_energy()
        if self._excitation_env is None:
            from .pyitensor import idmrg_excitations
            self._excitation_env = idmrg_excitations.build_excitation_environment(
                self._vumps_result)
        return self._excitation_env

    def excitation_energies(self, k, n=1):
        """The lowest `n` excitation energies (above the ground state) at
        momentum `k` (radians, per unit cell) of the tangent-space/
        quasiparticle excitation ansatz.

        itensor_version="python": see pyitensor.idmrg_excitations' own
        module docstring for the algorithm (requires gs_method="vumps",
        see `_get_excitation_environment`'s own docstring).

        itensor_version=3: C++ port of the same algorithm
        (mpscpp3/chain_session.h's Chain::vumps_excitation_energies) --
        also requires gs_method="vumps"; calls gs_energy() itself if not
        already run with that gs_method. See that C++ method's own doc
        comment for validation status (cross-checked directly against this
        same "python" path across a momentum scan on TFIM/Heisenberg
        models at D=1..3, matching to ~1e-10 or tighter -- see
        tests/test_vumps_excitations_v3.py).

        Any converged bond dimension D>=1 is supported on both backends."""
        if self.itensor_version == "python":
            env = self._get_excitation_environment()
            from .pyitensor import idmrg_excitations
            return idmrg_excitations.excitation_energies(env, k, n=n)
        if self.itensor_version != 3 or self.gs_method != "vumps":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.excitation_energies/excitation_gap: "
                "only itensor_version=\"python\" or itensor_version=3, both "
                "with gs_method=\"vumps\", are supported -- set "
                "self.gs_method = \"vumps\" before calling "
                "gs_energy()/excitation_energies()")
        if self._session3 is None or not self._session3_has_vumps:
            self.gs_energy()
        return np.array(self._session3.vumps_excitation_energies(k, n))

    def excitation_gap(self, ks=None):
        """The scalar "first excited state" gap: min_k
        excitation_energies(k, n=1) over a momentum scan (default
        `ks=numpy.linspace(-pi, pi, 41)`) -- mirrors the finite-chain
        `get_gap()` naming (docs/user_guide.md Sec.4), which returns
        `get_excited(n=2)[1]-get_excited(n=2)[0]`; here there is no
        discrete second state, only a momentum-resolved band, so the gap
        is defined as the band's own minimum. See excitation_energies'
        own docstring for the gs_method="vumps" requirement."""
        if ks is None:
            ks = np.linspace(-np.pi, np.pi, 41)
        return min(self.excitation_energies(k, n=1)[0] for k in ks)

    def _spectral_operator_matrix(self, opname, p):
        """The grouped-supersite d_g x d_g matrix for `opname` on sub-site
        `p` of the unit cell, in the M[i,o] convention every pyitensor
        infinite-chain contraction uses -- shared argument checking for
        spectral_weights/dynamical_structure_factor."""
        if not (0 <= p < self.n_uc):
            raise ValueError("spectral_weights: p must be in 0..{} (n_uc-1), "
                              "got {!r}".format(self.n_uc - 1, p))
        if is_fermionic(opname):
            # A parity-odd operator drags a Jordan-Wigner string that this
            # single-supersite ansatz would have to close at infinity -- the
            # same reason correlator() rejects an odd-parity operator pair.
            raise ValueError(
                "spectral_weights: {!r} is a fermionic (parity-odd) operator, "
                "whose Jordan-Wigner string cannot be closed inside a single "
                "unit cell -- only parity-even (e.g. density/spin) operators "
                "have a well-defined single-site spectral weight here".format(
                    opname))
        from .pyitensor import vumps
        result = self._vumps_result
        return vumps._embed_group_operator(
            result.sites_uc, result.n_uc,
            {p: result.sites_uc.site_type(p + 1).matrix(opname)})

    def spectral_weights(self, opname, k, p=0, n=1, return_total=False):
        """(energies, weights): the lowest `n` excitation energies at
        momentum `k` (radians per *unit cell*) together with the exact
        delta-peak spectral weight each of them carries for the local
        operator `opname` on sub-site `p` of the unit cell.

        `weights[a] = |<k,a| O(k) |Psi>|^2` with `O(k) = N^{-1/2} sum_m
        e^{ikm} O_m` and `|k,a>` the normalized quasiparticle state, i.e.
        the residue of the delta peak branch `a` contributes to the
        zero-temperature dynamical structure factor

            S(k, w) = sum_a weights[a] * delta(w - energies[a]).

        This is a genuine thermodynamic-limit, momentum-resolved dynamical
        correlator: unlike `kpm_finite`/`td_dynamical_correlator`, which
        embed a finite window in the infinite chain and therefore carry a
        window-size and a broadening/time-truncation error, the
        quasiparticle branches here are exact delta functions at exactly
        defined momenta. What they do *not* capture is multi-particle
        continuum weight, which is precisely what `return_total` measures
        (see below), so the two families of methods are complementary
        rather than one superseding the other.

        `return_total=True` appends a third return value: the total weight
        summed over EVERY branch, not just the `n` requested -- which is
        exactly the per-site connected static structure factor
        `sum_r e^{ikr} (<O_0 O_r> - <O>^2)` at this momentum, since a
        one-site operator applied to a uniform MPS lands exactly inside the
        variational tangent space (see pyitensor.idmrg_excitations.
        _spectral_source_vector's own "Sum rule" section). `weights.sum() /
        total` is then the fraction of the k-resolved spectral weight the
        returned branches carry -- the practical test of whether a
        single-mode picture is adequate at that momentum. It is free: the
        same source vector every weight is an inner product against.

        Backend support mirrors `excitation_energies`' own
        itensor_version="python" + gs_method="vumps" requirement, with one
        further restriction: itensor_version=3 is NOT supported here at all
        (excitation_energies does support it). The C++ port
        (Chain::vumps_excitation_energies) has neither the eigenvectors nor
        the mixed-transfer source vector this needs, the same way it has
        not yet picked up the iterative eigensolver -- see
        docs/idmrg_improvement_plan.md.

        Two things to know before reading an individual weight. Within a
        degenerate multiplet the split between branches is
        **basis-arbitrary** (the eigensolver picks an arbitrary basis of a
        degenerate eigenspace), so only multiplet sums are physical -- on
        the AKLT chain at `D=2` the branches split into an SU(2) triplet
        and a quintuplet at every momentum, and the triplet's three
        weights come out ~0.08/0.84/0.08 with their sum equal to the whole
        total. Away from an exactly-representable ground state the
        multiplet is only *approximately* degenerate -- split by the
        variational error rather than by machine precision (1.6e-5 on the
        Haldane chain at `D=16`, against 1e-13 at `D=12`, so not even
        monotonic in `D`) -- so group by the multiplet's size, not by a
        fixed energy tolerance. And a near-zero `weights.sum()/total` means the *wrong
        branches were asked for*, not that there is no response: on that
        same chain the two multiplets cross, so below `k~0.9` the lowest
        branch is the quintuplet, which `Sz` cannot reach at all (weight
        ~1e-23) while the triplet carrying every bit of the weight sits
        above it. Raise `n`; `return_total` is what makes this visible
        instead of silent.

        Sign convention: `k` labels the ansatz's own `sum_n e^{ikn}`; using
        the opposite sign relabels `k <-> -k` and changes nothing else."""
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.spectral_weights: only "
                "itensor_version=\"python\" is supported (the mpscpp3 port of "
                "the excitation ansatz returns energies only) -- got "
                "itensor_version={!r}".format(self.itensor_version))
        env = self._get_excitation_environment()
        from .pyitensor import idmrg_excitations
        M = self._spectral_operator_matrix(opname, p)
        return idmrg_excitations.spectral_weights(
            env, k, M, n=n, return_total=return_total)

    def dynamical_structure_factor(self, opname, ks=None, energies=None,
                                    delta=0.05, p=0, n=1):
        """(ks, energies, S): the quasiparticle contribution to the
        dynamical structure factor `S(k, w)` on a (momentum, energy) grid,
        Lorentzian-broadened by `delta`, built from `spectral_weights`'
        own exact delta peaks:

            S[i, j] = sum_a w_a(k_i) * (delta/pi)
                      / ((energies[j] - E_a(k_i))**2 + delta**2)

        `ks` defaults to `numpy.linspace(-pi, pi, 41)` (the same scan
        `excitation_gap` uses); `energies` defaults to 200 points spanning
        the computed band plus `5*delta` of margin on each side. `n` is how
        many branches per momentum to include -- with `n=1` this is the
        single-mode approximation to `S(k,w)`, exact for a model whose
        low-energy response is one isolated coherent mode and a *lower
        bound* on the true response otherwise (see `spectral_weights`'
        `return_total` for how to measure the shortfall).

        The broadening is purely cosmetic here, applied so the result can
        be plotted as a heat map -- the underlying peaks are exact delta
        functions at exact momenta, so `spectral_weights` is what to use
        for anything quantitative."""
        if ks is None:
            ks = np.linspace(-np.pi, np.pi, 41)
        ks = np.asarray(ks, dtype=float)
        bands, wgts = [], []
        for k in ks:
            e, w = self.spectral_weights(opname, k, p=p, n=n)
            bands.append(np.asarray(e, dtype=float))
            wgts.append(np.asarray(w, dtype=float))
        bands, wgts = np.array(bands), np.array(wgts)
        if energies is None:
            lo, hi = bands.min() - 5 * delta, bands.max() + 5 * delta
            energies = np.linspace(lo, hi, 200)
        energies = np.asarray(energies, dtype=float)
        gap = energies[None, :, None] - bands[:, None, :]
        lorentz = (delta / np.pi) / (gap ** 2 + delta ** 2)
        S = np.einsum('kea,ka->ke', lorentz, wgts)
        return ks, energies, S

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

        Requires `gs_method="idmrg"` on either backend. `window=0` is
        supported on both `itensor_version="python"` (pyitensor.idmrg.
        local_excitation_gap) and `itensor_version=3` (its C++ port,
        Chain::idmrg_local_excitation_gap). `window>0` is
        `itensor_version="python"` only: pyitensor.idmrg.
        local_excitation_gap_windowed is explicitly a prototype rather than
        stable API (see its own docstring), so it was deliberately left out
        of the v3 port rather than pinned down in C++ first."""
        if self.gs_method != "idmrg":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.local_excitation_gap: only "
                "gs_method=\"idmrg\" is supported -- see "
                "pyitensor.idmrg.local_excitation_gap's own docstring")
        if self.itensor_version == 3:
            if window != 0:
                raise NotImplementedError(
                    "Infinite_Many_Body_Chain.local_excitation_gap: "
                    "window>0 is supported only for "
                    "itensor_version=\"python\" -- "
                    "pyitensor.idmrg.local_excitation_gap_windowed is a "
                    "prototype, not stable API, and was deliberately not "
                    "ported to the v3 C++ backend")
            if self._session3 is None or self._session3_has_vumps:
                self.gs_energy()
            gap, e0_fresh, e0_stored = \
                self._session3.idmrg_local_excitation_gap_detail(niter)
            _warn_if_growth_missed_local_ground_state(e0_fresh, e0_stored)
            return gap
        if self.itensor_version != "python":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.local_excitation_gap: only "
                "itensor_version=\"python\" and itensor_version=3 are "
                "supported, got {!r}".format(self.itensor_version))
        if self._result is None:
            self.gs_energy()
        from .pyitensor import idmrg
        if window == 0:
            gap, e0_fresh, e0_stored = idmrg.local_excitation_gap(
                self._result, niter=niter, detail=True)
            _warn_if_growth_missed_local_ground_state(e0_fresh, e0_stored)
            return gap
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
        two_point_correlator's own r convention). Fermionic (parity-odd)
        names need no special handling here and get none: the finite window
        is an ordinary `Many_Body_Chain`, so its operators pass through
        `MultiOperator.to_terms()` and carry the *global*, finite-chain
        Jordan-Wigner string -- the correct convention for the finite object
        this method actually diagonalizes (pinned against exact free
        fermions in `tests/test_idmrg_window_fermionic.py`). This is
        unlike `td_dynamical_correlator`, which works with the infinite
        chain's own tensors directly and threads its own local string.

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
        method's own comment) are the two implemented backends; any other
        raises NotImplementedError.

        Both backends tile the *gauge-consistent* unit cell
        (`IDMRGResult.cell_raw` for `"python"`, `Chain::idmrg_cell_raw_`
        for v3), not the raw per-micro-step iDMRG factors: tiling the
        latter leaves the window's energy right while putting `S(x,t)`
        quantitatively wrong for *every* operator (it missed the exact
        `S(x,t=0) == correlator(...)` identity by up to ~1.7e-1 on a plain
        spin chain, with shape, finiteness and energy density all looking
        right, which is how it went unnoticed on v3 until 2026-08-29).
        Both now satisfy that identity to ~1e-10 or better. BOTH require
        `self.gs_method == "idmrg"`
        (raises NotImplementedError otherwise, mirroring
        `local_excitation_gap`'s own gate) -- unlike `vev`/`correlator`/
        `excitation_energies`, this feature has no `gs_method="vumps"`
        equivalent at all: it needs the growing algorithm's own converged
        environment snapshot (`IDMRGResult.env_HL`/`env_HR` for
        `"python"`, `Chain::idmrg_ground_state`'s internal snapshot for
        `itensor_version=3`), which VUMPS's mixed-gauge `{AL,AR,C,GL,GR}`
        has no equivalent of. Calls `gs_energy()` first if not already run
        under `gs_method="idmrg"` (`self._result`/`self._session3` unset
        or, for v3, populated by a stale `gs_method="vumps"` run instead --
        see `self._session3_has_vumps`). The v3 path's
        own `x_values` may not extend beyond the window's own explicit
        range (`center+x` must stay within `[1, n_window*n_uc]`) -- unlike
        the "python" backend, it does not pad beyond the window with extra
        unevolved unit-cell copies, so increase `n_window` instead if a
        wider `x_values` is needed.

        **Fermionic operators.** A parity-odd pair (`C`/`Cdag`, either
        order) computes the physical fermionic correlator: the Jordan-Wigner
        string is threaded on both sides -- explicitly across the window on
        the ket (before the evolution, so it is evolved along with the
        perturbation) and across the bra at measurement time -- exactly as
        `correlator` does for the static case, and pinned against it by the
        exact `S(x, t=0) == correlator(...)` identity in
        `tests/test_idmrg_window_fermionic.py`. See
        `pyitensor.idmrg_window.apply_local_operator` for where the two
        semi-infinite strings are truncated and why that truncation is the
        same boundary approximation the IBC window already makes. A pair
        with odd *total* fermion parity (one fermionic operator against a
        parity-even one) raises: its string can never close. `connected=True`
        subtracts nothing for a parity-odd pair, whose disconnected
        background `<A><B>` is zero by symmetry.

        `opname_j` is applied at sublattice position `p_i` (0..n_uc-1) and
        evolved forward in time; `opname_i` is inserted at the shifted
        position `p_i`'s own site `+x` (bra side, *not* time-evolved) --
        see `pyitensor.idmrg_window.dynamical_correlator_td`'s own
        docstring for the precise (Schrödinger-picture, background-
        subtracted, ground-state-energy-shifted) convention this follows,
        matching the rest of this codebase's own established "TD" submode
        convention (`mpscpp3::quench_tdvp`'s `Hshift = H - EGS*Id`).

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
            if self.gs_method != "idmrg":
                raise NotImplementedError(
                    "Infinite_Many_Body_Chain.td_dynamical_correlator: "
                    "itensor_version=\"python\" requires gs_method=\"idmrg\" "
                    "(needs pyitensor.idmrg.IDMRGResult's own env_HL/env_HR "
                    "snapshot, which gs_method=\"vumps\" does not build) -- "
                    "set self.gs_method = \"idmrg\" before calling "
                    "gs_energy()/td_dynamical_correlator()")
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
        if self.gs_method != "idmrg":
            raise NotImplementedError(
                "Infinite_Many_Body_Chain.td_dynamical_correlator: "
                "itensor_version=3 requires gs_method=\"idmrg\" (needs "
                "Chain::idmrg_ground_state's own converged-environment "
                "snapshot, which gs_method=\"vumps\" does not build) -- set "
                "self.gs_method = \"idmrg\" before calling "
                "gs_energy()/td_dynamical_correlator()")
        if self._session3 is None or self._session3_has_vumps:
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
        from .timedependent import sxt_to_skomega
        return sxt_to_skomega(ts, xs, S, dt, ks=ks, es=es, window=window,
                               delta=delta, factor=factor)


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
