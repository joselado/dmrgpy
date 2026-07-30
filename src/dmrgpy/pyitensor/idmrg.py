"""Infinite DMRG (iDMRG): the standard 2-site growing/infinite-size
algorithm (White, PRL 69, 2863 (1992)), generalized to an n_uc-site unit
cell per McCulloch's "Infinite size density matrix renormalization group,
revisited" (arXiv:0804.2509) -- following the shape of the reference
implementation at github.com/ITensor/iDMRG.

This module is self-contained within the pyitensor engine: it reuses
tensor.py/svd.py/kernels.py directly and imports only `_lanczos_ground_state`
from dmrg.py (a generic, chain-position-agnostic Lanczos solver), but never
imports `multioperator.py` or anything from the top-level dmrgpy package --
exactly mirroring how chain.py itself only ever consumes an already-
`MultiOperator.to_terms()`-shaped term list, never a MultiOperator object.
`infinitechain.py` (the Many_Body_Chain-level API) is responsible for
validating/canonicalizing the user's Hamiltonian into the plain
`MultiOperator.op`-shaped term lists (`h_intra_op`, `h_inter_op` below,
0-based site indices, no Jordan-Wigner threading needed -- this module
threads Jordan-Wigner strings itself, term by term, exactly like
autompo.HTerm.resolve() does, see `_term_site_matrices`) that this module's
functions take directly.

== Why the periodic MPO tensor is built directly, not extracted ==

An earlier version of this module built a small (m=3 unit cells) finite
reference system via the existing AutoMPO/mpobuilder.to_mpo machinery, and
tried to extract one repeating MPO tensor per sublattice position from its
middle repeat, relabeling bond Indices so the same tensor could be reused
indefinitely. Two independent bugs were found this way, both confirmed by
direct numerical checks (not just suspected):

1. `to_mpo`'s SVD-based compression picks an independent, arbitrary U/V
   basis at *every* bond, even at cutoff=0. Two positions of the compressed
   MPO that represent the exact same repeating physical bond can therefore
   come out expressed in different, non-interchangeable bases despite
   having the same dimension -- reusing one Index object across them (the
   whole point of the extraction) silently produced a measurably
   non-Hermitian effective 2-site Hamiltonian (~16% relative asymmetry)
   instead of raising. Building the reference MPO as an *exact,
   uncompressed* sum of per-term MPOs instead (no SVD at all) fixed this
   specific symptom -- but:
2. Even with that fixed, the reported energy density still converged to 0
   instead of the correct value. The root cause: a tensor extracted from
   any *fixed, finite* reference system has a bond dimension tied to that
   reference's own (finite) term count, with no channel that has an
   unconditional identity self-loop -- i.e. no genuine "accumulator" that
   can keep summing contributions across an *unbounded* number of repeats.
   Reusing such a tensor forever caps how much Hamiltonian content a fixed
   bond dimension can hold, which is exactly why the reported energy
   stopped growing extensively and instead decayed toward 0.

The fix implemented below is to build each sublattice's periodic MPO tensor
*directly*, as a genuine finite-state automaton (see `_build_automaton`):
one "F" (accumulator) channel with an Id self-loop, plus one "pending"
channel per still-open 2-site term -- the same structure any standard
finite-range-Hamiltonian MPO textbook derivation uses (e.g. the familiar
bond-dimension-5 NN Heisenberg MPO), generalized here to the "at most two
adjacent unit cells" term range `infinitechain.py` validates. Because this
tensor is built directly rather than extracted from an independently
gauge-fixed system, its own left and right bond Indices are, by
construction, indexed by the identical, explicitly-enumerated channel list
-- no relabeling, no gauge ambiguity, and (unlike the finite reference) an
unconditional identity self-loop on the accumulator channel, so the growth
loop can sum an unbounded number of repeats at fixed bond dimension.

== The growing algorithm, concretely ==

1. Build the periodic per-sublattice automaton tensors W[0..n_uc-1] (see
   `_build_automaton`, channel space "S" + "F" (accumulator) + pending),
   plus two dedicated *boundary* tensors, both via `_project_channel`:
   project W[0]'s left leg onto "S" for the very first site ever absorbed
   into HL, and W[n_uc-1]'s right leg onto "F" for the very first site
   ever absorbed into HR.
2. Grow two environments HL/HR, cold-started at `HL=HR=None` (the same
   sentinel dmrg.py's own `_all_left_environments`/`_all_right_environments`
   use for an open chain's true boundary). Each *macro-iteration* performs
   `n_uc` *micro-steps* (one full unit cell added to each side); at
   micro-step `m` (0..n_uc-1) the two newly-inserted sites have sublattice
   `p_L=m` (extending HL) and `p_R=n_uc-1-m` (extending HR). This indexing
   keeps each side's *own* absorption order correctly continuous (HL always
   absorbs 0,1,...,n_uc-1 in order; HR's own internal order, after
   prepending, works out the same way -- verified by hand for n_uc=2), but
   the two *active* sites phys_L/phys_R of one micro-step are only
   genuinely adjacent to each other (a real requirement for the 2-site
   solve below to be physically meaningful) when n_uc<=2 -- for n_uc>=3 an
   intermediate micro-step's phys_L and phys_R are several sites apart with
   real, not-yet-inserted sites in between, and `_build_automaton`'s own
   per-sublattice pending-channel content (see its docstring) differs
   between the two, so there is no way to even wire them together
   correctly without a genuine redesign of this pairing scheme. `n_uc>=3`
   is therefore rejected explicitly (`infinitechain.py`'s constructor),
   rather than left to fail confusingly deep inside a Lanczos call.
3. Each micro-step's 2-site local ground state is found by the same
   matvec-and-Lanczos machinery dmrg.py's own two_site_heff/
   _local_ground_state use (`kernels.make_matvec` + `dmrg._lanczos_ground_state`,
   imported directly, not duplicated), SVD-split (`svd.svd`) into a new
   left-canonical tensor (absorbed into HL) and right-canonical tensor
   (absorbed into HR). `idmrg_ground_state` re-mints a fresh mpo-axis
   Index for HL's and HR's own environment tensor at the end of *every*
   micro-step (see its own inline comment) rather than reusing whatever
   Index `_build_automaton` happened to label that channel space with --
   required because `_build_automaton`'s `boundary_idx` pool has only
   `n_uc` distinct Index objects total, reused verbatim every time a given
   sublattice position comes up again (every macro-iteration past the
   first); left un-broken, HL's and HR's own mpo axes can end up being the
   *same* Index object (guaranteed for n_uc=1, where a single sublattice's
   left and right automaton legs are literally identical) once both
   environments are non-trivial and used together in the same
   `_local_two_site_solve` call, corrupting the effective Hamiltonian from
   the second macro-iteration onward -- confirmed by extracting the dense
   Heff at macro-iteration 2 of a uniform n_uc=1 chain and finding its full
   eigenvalue spectrum did not match exact 6-site ED at all (not a
   truncation effect: macro-iteration 1's spectrum matched ED to machine
   precision).
4. Convergence: once per macro-iteration, the finite-difference of the
   local ground-state eigenvalue between consecutive macro-iterations,
   divided by 2*n_uc (exactly the number of new physical sites per
   macro-iteration), is compared against `etol` -- the standard infinite-
   algorithm energy-density diagnostic. This is an *energy* criterion only:
   see `_local_two_site_solve`'s and `IDMRGResult.state_overlap`'s own
   docstrings for why energy-density convergence alone does not guarantee
   `U_list` itself has settled into a self-consistent, translationally-
   invariant unit cell (a real, previously-unfixed bug for gapless/
   SU(2)-symmetric models -- fixed by warm-starting each macro-iteration's
   local solve from the previous one's own local ground vector, see point
   3's own solve step and `state_overlap` for the diagnostic that confirms
   it).

== Static correlators ==

After convergence, the last macro-iteration's own per-micro-step
left-canonical tensors (one per sublattice position) are taken as the
converged unit cell's canonical MPS tensors, and one-/two-point static
expectation values are computed via the standard infinite-MPS
transfer-matrix formalism (dominant left/right fixed points of the
per-unit-cell transfer operator) -- see `onsite_expectation`/
`two_point_correlator`.
"""

import numpy as np

from . import kernels
from .dmrg import _lanczos_ground_state
from .index import Index
from .sites import SiteX
from .sites.base import is_fermionic
from .svd import _truncate, svd
from .tensor import ITensor, contract_many, dag
from .tensor import noPrime as _t_noPrime
from .tensor import prime as _t_prime


def _unprimed_site_index(W):
    """The unprimed ("ket"-side) Site-tagged Index of an MPO tensor W --
    always its own physical leg, regardless of whether W is a bulk tensor
    or one of the two projected boundary tensors (see `_project_channel`),
    which is why this is read off W itself rather than looked up
    separately."""
    return next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 0)


def _fresh_physical_copy(W):
    """W with a freshly-minted (unprimed,primed) physical Index pair,
    every other (Link/mpo-bond) leg untouched -- needed whenever the exact
    same per-sublattice MPO tensor object is used for *both* active sites
    of one micro-step (p_L==p_R, which happens at exactly one micro-step
    per macro-iteration whenever n_uc is odd, including every micro-step
    of a uniform n_uc=1 chain). Without this, the "two different" physical
    legs of the local 2-site problem would collide onto the very same
    Index object (since they'd both be read off the identical W tensor),
    silently corrupting the local eigenproblem instead of raising -- caught
    directly via a ValueError from kernels.py's index-matching machinery
    on an n_uc=1 test case before this fix was added."""
    old = _unprimed_site_index(W)
    new = old.sim()
    new_inds = tuple(new.setprime(ind.plev) if ind.id == old.id else ind
                      for ind in W.inds)
    return ITensor(new_inds, W.array)


def _term_site_matrices(op_term, sites_uc, n_uc):
    """op_term: (coef, [(name,rel_site), ...]) in MultiOperator.op format,
    rel_site 0-based, spanning at most 2 distinct sites (validated by
    infinitechain.py's `set_hamiltonian`). Returns (rel_sites_sorted, coef,
    mats, ferm): rel_sites_sorted is the sorted list of distinct touched
    sites (length 1 or 2); mats[k] is the composed (dim,dim) matrix at
    rel_sites_sorted[k] (multiple operator factors at the same site are
    combined via matrix product, reversed and un-transposed relative to
    the term's own given order -- the same net composition
    autompo.HTerm.resolve() computes for its "stored", (in,out)-convention
    per-site matrix, worked out by tracing through its own std-convention
    `.T`/composition steps); ferm[k] is whether an odd number of fermionic
    factors occur at that site -- used by `_classify_terms` to decide
    whether the "still open" channel between two touched sites of a 2-site
    term must propagate a Jordan-Wigner "F" string or a plain identity
    (mirrors HTerm.resolve()'s own `is_site_fermionic`-triggered carry
    flag, which only ever depends on the *first* touched site of a term)."""
    coef = op_term[0]
    by_site = {}
    for name, site in op_term[1:]:
        by_site.setdefault(site, []).append(name)
    rel_sites = sorted(by_site)
    mats, ferm = [], []
    for site in rel_sites:
        names = by_site[site]
        st = sites_uc.site_type((site % n_uc) + 1)
        combined = st.matrix(names[-1])
        for nm in reversed(names[:-1]):
            combined = combined @ st.matrix(nm)
        mats.append(combined)
        ferm.append(sum(1 for nm in names if is_fermionic(nm)) % 2 == 1)
    return rel_sites, complex(coef), mats, ferm


def _classify_terms(h_intra_op, h_inter_op, sites_uc, n_uc):
    """h_intra_op/h_inter_op -> (onsite, bonds): onsite is a list of
    (rel_site, coef, mat) 1-site terms; bonds is a list of dicts
    {rel_a, mat_a, rel_b, mat_b, carry_ferm, coef} for each 2-site term
    (rel_a<rel_b, rel_a always < n_uc -- guaranteed by
    infinitechain.py's canonicalization, which shifts any "pure R" term,
    touching only site indices >= n_uc, back down into h_intra)."""
    onsite, bonds = [], []
    for op_term in list(h_intra_op) + list(h_inter_op):
        rel_sites, coef, mats, ferm = _term_site_matrices(op_term, sites_uc, n_uc)
        if len(rel_sites) == 1:
            onsite.append((rel_sites[0], coef, mats[0]))
        elif len(rel_sites) == 2:
            a, b = rel_sites
            if a >= n_uc:
                raise ValueError(
                    "idmrg: internal error -- inter-cell term does not "
                    "touch the central cell ({})".format(op_term))
            bonds.append(dict(rel_a=a, mat_a=mats[0], rel_b=b, mat_b=mats[1],
                               carry_ferm=ferm[0], coef=coef))
        else:
            raise ValueError(
                "idmrg: terms touching more than 2 distinct sites are not "
                "supported ({})".format(op_term))
    return onsite, bonds


def _active_channels_at(bonds, n_uc, p):
    """[(bond_index, r), ...] -- every 2-site term with a "pending"
    instance active just before absorbing sublattice p (i.e. this is the
    channel list valid *entering* sublattice p's own tensor, having
    already absorbed everything up to but not including it), where r
    counts down the number of *more* sites needed (including the one
    about to be absorbed) until that specific instance completes: r=reach
    means "just absorbed its own rel_a, about to absorb rel_a+1" (the
    first position where it's actually pending -- rel_a's own site itself
    is *not* pending yet when entering it, only starting in its output);
    r=1 means "completes upon absorbing this very site". A term instance
    with window [rel_a, rel_b] is pending entering absolute position P
    (P mod n_uc = p) iff rel_a < P <= rel_b (mod the unit-cell tiling),
    with r = rel_b - P + 1 -- solving for which p a given r corresponds
    to gives p = (rel_b - r + 1) % n_uc, the formula used below (an
    earlier version omitted the "+1", silently checking a mismatched
    sublattice for every n_uc > 1 -- invisible for n_uc=1, where
    (p+1)%1 == p makes the two formulas coincide, but fatal for n_uc=2:
    it produced a completely decoupled, exactly-zero effective
    Hamiltonian instead of raising). A bond can appear more than once
    here (distinct r values) if its reach exceeds n_uc, meaning two
    different unit cells' copies of the same term have genuinely
    overlapping pending windows -- within this module's documented scope
    (couplings spanning at most two adjacent unit cells) this happens for
    at most 2 simultaneous instances of any one term."""
    out = []
    for bi, b in enumerate(bonds):
        reach = b["rel_b"] - b["rel_a"]
        for r in range(1, reach + 1):
            if (b["rel_b"] - r + 1) % n_uc == p:
                out.append((bi, r))
    return out


def _onsite_matrix(onsite_by_p, p, d):
    """Σ coef*mat over every onsite term attributed to sublattice p, or
    None if there is none (kept separate from np.zeros so callers can
    tell "no onsite term here" from "an onsite term that happens to be
    the zero matrix", and so the "F->F" self-loop's mandatory Id doesn't
    need a redundant +0).

    `onsite_by_p[p]` holds (coef, mat) pairs -- `rel_site` (`_classify_terms`'
    own 3rd list entry) is already consumed as this dict's own key by
    `_build_periodic_mpo`, not carried into each stored tuple. A previous
    version of this loop unpacked 3 values here (`for _rel, coef, m in ...`)
    -- a leftover from `_classify_terms`' own (rel_site, coef, mat) shape
    that was never updated to match -- which raised `ValueError: not enough
    values to unpack` for *any* Hamiltonian containing an onsite term (bond-
    only Hamiltonians never populate onsite_by_p, so every existing test
    happened to avoid this path); confirmed directly and fixed here."""
    mat = None
    for coef, m in onsite_by_p.get(p, []):
        mat = coef * m if mat is None else mat + coef * m
    return mat


def _build_automaton(h_intra_op, h_inter_op, site_types, n_uc):
    """sites_uc = SiteX(site_types), then _build_periodic_mpo(...) against
    it -- split out so apply_mpo's own callers can build a *different*
    (bounded, non-Hamiltonian) automaton sharing an *existing* IDMRGResult's
    own sites_uc (required for the physical Indices to match by identity,
    see _build_periodic_mpo's own docstring), rather than always minting a
    fresh SiteX the way a Hamiltonian build always wants to."""
    sites_uc = SiteX(list(site_types))
    return sites_uc, _build_periodic_mpo(h_intra_op, h_inter_op, sites_uc, n_uc)


def _build_periodic_mpo(h_intra_op, h_inter_op, sites_uc, n_uc):
    """Build the periodic per-sublattice MPO tensors W_bulk[0..n_uc-1]
    directly, as a genuine finite-state automaton -- see this module's
    docstring for why this, not an extraction from any fixed finite
    reference system, is required.

    Two *distinct* trivial channels, both persistent (part of every
    W_bulk[p]'s own reused channel space, not a one-off boundary
    construct), plus one "pending" channel per still-open 2-site term:

    - "S" (index 0): self-loops via Id, and is the *only* channel a new
      2-site term may start from. Nothing ever transitions into S from
      anywhere else, so it always carries pure, unweighted Id content --
      never contaminated by whatever has already been accumulated. Also
      the *only* channel a site's own onsite term may start from, via a
      direct S->F transition (see the "S,F" case below) -- an onsite term
      is, in automaton terms, a "bond" that starts and completes on the
      very same site.
    - "F" (index 1, the accumulator): self-loops via plain Id (nothing
      added -- see below for why), and receives completed pending
      channels *and* direct onsite completions from S. Deliberately
      *cannot* start a new term itself.

    A site's own onsite term must live on the direct "S,F" transition, NOT
    added onto "F,F"'s own self-loop (an earlier, wrong version of this
    function did exactly that -- `mat = Id; if onsite_mat is not None: mat
    += onsite_mat` on the F,F entry). Two independent things go wrong with
    putting it there instead of on S,F:

    1. No branch below ever actually sets the "S,F" entry (every existing
       branch requires `rch`/`lch` to be a pending-channel tuple, never
       bare "S" or "F"), so with F,F holding the onsite content, W stays
       *block-diagonal* between {S} and {F,pending...} for any Hamiltonian
       with no bond terms at all (confirmed directly: S,F and F,S both
       stay exactly 0 in that case) -- and a matrix power/product of a
       block-diagonal matrix stays block-diagonal forever, so <S|W^N|F>
       (exactly what the growing algorithm's boundary projections extract,
       see `_project_channel`) is identically 0 for every N, silently
       dropping every onsite term entirely. Confirmed directly: a bare
       field Hamiltonian (`-B*SzC[0]`, n_uc=1, no bond term) converged to
       energy density 0 and <Sz>~0 regardless of B, instead of the exact
       B/2 (fully polarized, since a field-only Hamiltonian is exactly
       solvable).
    2. Once at least one bond term *does* eventually activate F (a
       completed pending channel), every *further* site absorbed while F
       is already hot picks up "Id+onsite_mat" again on top of whatever HL
       already represents -- but HL's own recursive `_extend_HL` already
       carries the *entire* accumulated-so-far MPO-chain product as a live
       tensor leg (not a collapsed scalar), so this re-adds that site's
       onsite content into an already-summed F channel instead of summing
       it in exactly once via its own site's S->F transition -- an
       exponential (multiplicative-in-effect) blow-up, not a linear
       (additive) one. Confirmed directly: an n_uc=2 XXZ chain with a
       small field (`-0.1*(SzC[0]+SzC[1])` added to an otherwise-normal
       bond Hamiltonian) reported energy densities of -1.3e7, then -1.4e23
       at B=0.3, -3.6e37 at B=0.5, -1.3e69 at B=1.0 -- diverging, not
       converging, with `.converged` staying False throughout -- instead
       of a finite, B-dependent shift of the B=0 energy density. The fix
       (S,F direct transition carrying onsite_mat, plain Id with no
       addition on F,F) is exactly the standard textbook automaton-MPO
       construction for summing an onsite term into a running total (the
       same upper-triangular-in-channel-index structure any finite-range
       Hamiltonian MPO derivation uses, e.g. the familiar bond-dimension-5
       NN Heisenberg-plus-field MPO's own field entry) -- verified by hand
       for the matrix-product identity this construction relies on
       (an upper-triangular [[Id,h],[0,Id]]-per-site chain product's own
       top-right entry sums exactly sum_i h_i, by induction on chain
       length) and confirmed numerically: the same field test above now
       reproduces the field=0 energy density plus the expected
       field-induced shift, converging cleanly for every B tried.

    Two earlier, both wrong, designs were tried and are worth recording
    here since the failure modes were not obvious in advance:

    1. A single merged S/F channel (this module's very first version):
       the "nothing has ever happened, anywhere" path then trivially
       reaches the same index used to terminate the chain, contributing a
       spurious, wavefunction-coupled identity-operator term -- confirmed
       directly as a uniform +1 shift of the exact 2-site spectrum for a
       simple test case.
    2. Separate S and F, but S consumed entirely by a one-off boundary
       tensor (not part of W_bulk's own reused structure), *and* F itself
       allowed to also start new pending channels (mirroring S's own
       "start" edge, reasoning that F should be able to start a *second*
       independent bond once a first one has already completed). This is
       wrong for a different reason: starting a new term from F multiplies
       the new term by *whatever F already holds* (a product), not an
       independent, additive contribution -- the correct way to add a new,
       independent term to a running sum is for it to start from something
       carrying pure Id weight, i.e. S specifically, not F. Confirmed
       directly: this version could not represent a growth step's *newest*
       bond at all whenever F was still zero (e.g. every micro-step of the
       very first macro-iteration, where nothing has completed yet), and
       once S was instead made persistent (to fix that), letting *both* S
       and F start new terms let every later bond be double-counted
       through two structurally distinct paths that never resolve against
       each other, since the growth loop has no final right boundary to
       exclude the "still S" branch -- reported energy diverged to
       -infinity across iterations instead of converging.

    The fix implemented below: S persists (self-loop, can start new
    terms), F persists (self-loop with onsite pickup, receives
    completions) -- but *only* S may start a new pending channel.

    `sites_uc` (an already-built SiteX -- see `_build_automaton`, the
    Hamiltonian-facing wrapper that mints a fresh one; apply_mpo's own
    bounded-operator callers instead pass an *existing* IDMRGResult's own
    sites_uc, required so the physical Indices below match, by identity,
    the converged U_list's own physical legs -- see idmrg.py's
    apply_mpo/grow_by_mpo docstrings) is where physical-leg Indices and
    per-site dimensions are looked up from, and h_intra_op/h_inter_op are
    resolved against its own named-operator matrices.

    Returns W_bulk: W_bulk[p] (rank 4: left-bond, right-bond, phys, phys')
    is sublattice p's tensor, with W_bulk[p]'s right bond and
    W_bulk[(p+1)%n_uc]'s left bond always the *same* canonical Index
    object -- valid by construction here (both are indexed by the
    identical, explicitly-enumerated channel list; there is no independent
    gauge choice anywhere in this construction to make them inconsistent,
    unlike extracting from an SVD-compressed reference)."""
    onsite, bonds = _classify_terms(h_intra_op, h_inter_op, sites_uc, n_uc)
    onsite_by_p = {}
    for rel, coef, mat in onsite:
        onsite_by_p.setdefault(rel, []).append((coef, mat))

    # chans[p][0] = "S", chans[p][1] = "F", chans[p][k>1] = a
    # (bond_index, r) pending-channel tuple -- the full channel list valid
    # just before absorbing sublattice p, shared as *both* W_bulk[p]'s own
    # left bond and W_bulk[(p-1)%n_uc]'s own right bond.
    S, F = "S", "F"
    chans = [[S, F] + _active_channels_at(bonds, n_uc, p) for p in range(n_uc)]
    boundary_idx = [Index(len(chans[p]), tags="Link") for p in range(n_uc)]
    dims = [sites_uc.dim(p + 1) for p in range(n_uc)]

    W_bulk = []
    for p in range(n_uc):
        left, right = chans[p], chans[(p + 1) % n_uc]
        left_idx, right_idx = boundary_idx[p], boundary_idx[(p + 1) % n_uc]
        d = dims[p]
        st = sites_uc.site_type(p + 1)
        s = sites_uc.si(p + 1)
        onsite_mat = _onsite_matrix(onsite_by_p, p, d)

        arr = np.zeros((len(left), len(right), d, d), dtype=complex)
        for li, lch in enumerate(left):
            for ri, rch in enumerate(right):
                mat = None
                if lch == S and rch == S:
                    mat = np.eye(d, dtype=complex)
                elif lch == F and rch == F:
                    mat = np.eye(d, dtype=complex)
                elif lch == S and rch == F:
                    # this site's own onsite term: starts and completes in
                    # one step, direct from S into the accumulator (see
                    # this function's own docstring for why this, not
                    # adding onsite_mat onto the F,F self-loop, is correct).
                    if onsite_mat is not None:
                        mat = onsite_mat
                elif lch == S and rch not in (S, F):
                    # a new pending channel starts here -- from S only
                    # (never from F, see this function's own docstring)
                    # -- iff this is exactly its own rel_a.
                    bi, r = rch
                    b = bonds[bi]
                    if r == b["rel_b"] - b["rel_a"] and b["rel_a"] == p:
                        mat = b["coef"] * b["mat_a"]
                elif rch == F and lch not in (S, F):
                    # a pending channel completes into the accumulator.
                    bi, r = lch
                    if r == 1:
                        mat = bonds[bi]["mat_b"]
                elif lch not in (S, F) and rch not in (S, F):
                    # a pending channel propagates one more site (Id, or
                    # "F" if this term's own Jordan-Wigner string is open).
                    bi_l, r_l = lch
                    bi_r, r_r = rch
                    if bi_l == bi_r and r_r == r_l - 1:
                        b = bonds[bi_l]
                        mat = st.matrix("F") if b["carry_ferm"] else np.eye(d, dtype=complex)
                if mat is not None:
                    arr[li, ri] = mat
        T = ITensor((left_idx, s, s.prime(1), right_idx),
                     np.transpose(arr, (0, 2, 3, 1)))
        W_bulk.append(T)
    return W_bulk


def _project_channel(T, side, idx):
    """T with its 'side' ('left' or 'right') bond leg contracted against
    the unit vector e_idx -- gives the true open-chain boundary tensor
    (rank 3) needed to cold-start the growth loop's very first micro-step:
    idx=0 ("S") for the left boundary (the chain has nothing accumulated
    yet), idx=1 ("F") for the right boundary (the chain reports only its
    genuinely-accumulated total, excluding the -- structurally impossible
    to reach from S alone, see `_build_automaton`'s docstring -- "S was
    still never resolved" case), exactly mirroring how a genuine finite
    chain's own leftmost/rightmost MPO tensor has no dangling boundary leg
    at all.

    Implemented via direct positional array indexing (np.take), not a
    generic Index-identity-matched ITensor contraction (`T * e`): for
    n_uc=1 (and more generally whenever a bulk tensor's left and right
    bond legs happen to be the very same canonical Index object, which
    `_build_automaton` deliberately arranges whenever a channel type
    recurs with period n_uc -- see its docstring), `T`'s left and right
    legs are indistinguishable *by identity* to the generic `mul_plan`
    matcher, which greedily matches whichever occurrence it scans first
    (always the left leg, since it's T.inds[0]) regardless of which side
    was actually requested -- confirmed directly: projecting side="right"
    for n_uc=1 silently sliced the *left* axis instead, giving a tensor
    representing "start a new term here" on both boundaries at once
    instead of one "start" and one "end", which produced an identically
    zero effective Hamiltonian for the very first growth step. Selecting
    the axis by its *position* in T.inds (side="left" -> axis 0,
    side="right" -> axis -1) rather than by Index identity sidesteps this
    ambiguity entirely."""
    axis = 0 if side == "left" else -1
    new_inds = T.inds[1:] if side == "left" else T.inds[:-1]
    new_array = np.take(T.array, idx, axis=axis)
    return ITensor(new_inds, new_array)


def _relabel_pos(T, pos, new_ind):
    """T with the Index at axis position `pos` replaced by `new_ind`,
    every other axis (including the array data) untouched. Selected by
    *position*, not by matching the old Index's identity (the same
    reasoning as `_project_channel`'s own docstring): whenever the axis
    being replaced shares its Index object with another axis of the same
    tensor -- always true of W_bulk[p]'s own left/right legs for n_uc=1,
    and possible for other small n_uc -- an identity-based relabel would
    silently rewrite *both* occurrences instead of just the one at `pos`.
    Used by `idmrg_ground_state` to give HL's and HR's own mpo axes an
    identity independent of `_build_automaton`'s small, reused
    boundary_idx pool -- see its own comment for why that reuse is
    otherwise unsafe."""
    inds = list(T.inds)
    inds[pos] = new_ind
    return ITensor(tuple(inds), T.array)


def _extend_HL(HL, HL_bra, W_p, U, left_ket_old, right_ket_new):
    """Absorb the newly-solved left-canonical site tensor U into HL, using
    MPO tensor W_p for the sublattice just absorbed. `left_ket_old` is U's
    own link toward the *existing* HL (None only for the very first-ever
    absorption); `right_ket_new` is U's own freshly-minted link away from
    HL (always present -- this becomes the new HL's own ket-side leg).
    Mirrors dmrg.py's `_extend_left`/`_relabel_bra_local`, but expressed
    directly in terms of already-known link Indices (this module has no
    `_Chain` container to search via `_link_at`, since the growing
    environment isn't attached to any fixed-length chain object) -- same
    contract_many()-routed, bra/ket-link-relabeled construction either way.
    Returns (new_HL, new_HL_bra) where new_HL_bra is right_ket_new's own
    bra-side counterpart, to be threaded into the *next* call."""
    right_bra_new = right_ket_new.sim()
    new_inds = []
    for ind in U.inds:
        if left_ket_old is not None and ind == left_ket_old:
            new_inds.append(HL_bra)
        elif ind == right_ket_new:
            new_inds.append(right_bra_new)
        else:
            new_inds.append(ind)
    bra_piece = dag(_t_prime(ITensor(tuple(new_inds), U.array), "Site"))
    pieces = [p for p in (HL, bra_piece, W_p, U) if p is not None]
    new_HL = contract_many(pieces)
    return new_HL, right_bra_new


def _extend_HR(HR, HR_bra, W_p, V, right_ket_old, left_ket_new):
    """Mirror of `_extend_HL`: V (the newly-solved right-canonical site
    tensor) is absorbed into HR by prepending it on HR's left. Returns
    (new_HR, new_HR_bra)."""
    left_bra_new = left_ket_new.sim()
    new_inds = []
    for ind in V.inds:
        if right_ket_old is not None and ind == right_ket_old:
            new_inds.append(HR_bra)
        elif ind == left_ket_new:
            new_inds.append(left_bra_new)
        else:
            new_inds.append(ind)
    bra_piece = dag(_t_prime(ITensor(tuple(new_inds), V.array), "Site"))
    pieces = [p for p in (HR, bra_piece, W_p, V) if p is not None]
    new_HR = contract_many(pieces)
    return new_HR, left_bra_new


def _local_two_site_solve(HL, HL_bra, HL_ket, W_pL, phys_L,
                           W_pR, phys_R, HR, HR_bra, HR_ket,
                           cutoff, maxdim, niter, x0_warm=None):
    """One micro-step's local ground-state solve: the effective 2-site
    Hamiltonian sandwiched by (HL, W_pL, W_pR, HR), diagonalized via the
    same matvec+Lanczos machinery as dmrg.py's own two_site_heff/
    _local_ground_state (see this module's docstring). Unlike finite DMRG's
    ket-tensor warm start (`two_site_heff`'s own `x0`), there is no
    persistent ket object here to read a warm start off directly -- every
    micro-step inserts brand-new *physical sites*, so there is no "this
    exact tensor, one sweep ago" to reuse. But the *shape* of the local
    eigenproblem at a given unit-cell position is identical macro-iteration
    to macro-iteration once bond dimension saturates at `maxdim`, so the
    caller (`idmrg_ground_state`) instead threads through the previous
    macro-iteration's own flattened ground vector at this same position as
    `x0_warm`, reused here as the Lanczos start whenever its size matches
    the current local dimension (falling back to a fresh random vector
    otherwise, e.g. before bond dimension has saturated). This is not
    cosmetic: an earlier, always-fresh-random version of this function let
    Lanczos land on an arbitrary member of a (near-)degenerate local ground
    manifold every macro-iteration (routine for gapless/SU(2)-symmetric
    models) -- the reported *energy* still converged fine (a degenerate
    manifold shares one eigenvalue) but `U_list` kept jumping between
    different members of that manifold instead of settling into one
    self-consistent, translationally-invariant state, corrupting every
    downstream static correlator (confirmed directly: the <H_uc>
    self-consistency identity, which any genuinely converged unit cell must
    satisfy exactly, was off by 0.4-0.6 -- not shrinking with `maxiter`, and
    not fixed by best-of-6 independent random-seed reruns, ruling out
    ordinary seed noise as the explanation). Warm-starting instead biases
    each macro-iteration's solve toward continuity with the previous one,
    letting the sequence of local ground states converge to a fixed member
    of the manifold. Returns (energy, U, S, V, evec0) -- the SVD split of
    the local ground state (truncated to (cutoff, maxdim)) plus the raw,
    un-truncated flattened ground vector for the *next* macro-iteration's
    own warm start."""
    order_in = ([HL_ket] if HL_ket is not None else []) + [phys_L, phys_R] + \
               ([HR_ket] if HR_ket is not None else [])
    shape = tuple(ind.dim for ind in order_in)
    s_L_out, s_R_out = phys_L.prime(1), phys_R.prime(1)
    order_out = ([HL_bra] if HL_bra is not None else []) + [s_L_out, s_R_out] + \
                ([HR_bra] if HR_bra is not None else [])
    pieces = [p for p in (HL, W_pL, W_pR, HR) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    dim = int(np.prod(shape))
    if dim <= 3:
        # too small a space for Lanczos to be meaningful; diagonalize directly
        # (mirrors dmrg.py's own _local_ground_state small-dim fallback).
        basis = np.eye(dim, dtype=complex)
        Hmat = np.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, v = np.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        eval0, evec0 = w[0], v[:, 0]
    else:
        # x0_warm is a flat array positionally aligned with `order_in`
        # (see idmrg_ground_state's own comment) -- its size, not the
        # identity of the Index objects behind `shape`, is what matters
        # here (the same positional-reuse convention this module's
        # _project_channel/_relabel_pos already rely on for other reasons).
        if x0_warm is not None and x0_warm.size == dim:
            norm = np.linalg.norm(x0_warm)
            v0 = x0_warm / norm if norm > 0 else None
        else:
            v0 = None
        if v0 is None:
            rng = np.random.default_rng()
            v0 = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        # Same niter floor as dmrg.py's _local_ground_state -- flagged there
        # as load-bearing (dropping it caused real, measured convergence
        # regressions); still applied here even with a warm start, since a
        # warm start only biases *which* near-degenerate solution Lanczos
        # converges to, not how many iterations that convergence itself
        # needs.
        eval0, evec0 = _lanczos_ground_state(matvec, v0, niter=max(niter, 200))

    theta = ITensor(tuple(order_in), evec0.reshape(shape))
    left_inds = ([HL_ket] if HL_ket is not None else []) + [phys_L]
    U, S, V, spec = svd(theta, left_inds, cutoff=cutoff, maxdim=maxdim)
    # U and V are used *without* S: HL/HR are block operators built by a
    # similarity transform through U (or V) alone -- <dag(U)|W_p|U>, in
    # _extend_HL/_extend_HR -- exactly mirroring White's original
    # infinite-algorithm block update (and dmrg.py's own _extend_left/
    # _extend_right, which likewise only ever see one SVD factor per
    # extension, via ket.A(i) after _apply_local_update has already
    # absorbed S into the *other* site). A version of this function that
    # instead absorbed sqrt(S) into both U and V (reasoning that this
    # growth loop extends HL and HR simultaneously, so neither should
    # asymmetrically "own" S) was tried and is wrong: multiplying a new
    # bond's own starting weight by whatever a neighboring environment's
    # bond matrix already contains, instead of the environment being a
    # clean similarity-transformed operator, corrupted every step past
    # the very first one -- confirmed directly, it broke the variational
    # monotonicity property (macro-iteration energies must only ever
    # decrease as more sites are absorbed) by the *second* macro-iteration
    # of a simple test case, and using plain U/V (this function's current
    # form) restores both exact agreement with independent ED at small
    # sizes and strict monotonicity at every iteration checked.
    return float(eval0.real), U, S, V, evec0


class IDMRGResult:
    """Converged iDMRG state, as needed for static correlators: the
    reference SiteX (for named-operator matrix lookups) and the last
    macro-iteration's own per-sublattice left-canonical site tensors
    (U_list[p], p=0..n_uc-1) -- a good approximation to the converged
    infinite chain's own canonical unit-cell MPS tensors, since by the
    last macro-iteration HL/HR are themselves converged to (near-)
    translational invariance.

    `state_overlap` is a diagnostic, not used by `.converged` (which
    remains a pure energy-density criterion, unchanged): the smallest
    (worst) per-sublattice normalized overlap
    `abs(<evec0_prev|evec0_new>)` between the final macro-iteration's own
    warm-started local ground vector (see `_local_two_site_solve`) and the
    one it was warm-started from, across all n_uc unit-cell positions.
    Close to 1 means the local ground states have actually stopped
    changing (a genuinely self-consistent, translationally-invariant unit
    cell -- the condition static correlators need), not just that the
    energy density has. `None` if unavailable (e.g. `niter_done` too
    small for every position to have been warm-started at least once).

    `local_superblock` is the raw ingredients of the very last micro-step's
    2-site local eigenproblem actually solved (HL/HR environments, the two
    MPO tensors, the two physical Indices, and the converged local ground
    vector+energy) -- kept around so `local_excitation_gap` can re-diagonalize
    that *same* effective Hamiltonian for its second-lowest eigenpair without
    having to rerun the growing algorithm. Not meant to be used directly by
    callers outside this module.

    `W_bulk`/`window_boundary` expose the converged *environment Hamiltonian*
    blocks -- the growth loop's own HL/HR as they stood entering the last
    executed macro-iteration (see idmrg_ground_state's own
    `env_window_boundary` comment for why exactly this snapshot, not the
    final post-loop HL/HR, is the one consistent with U_list's own edge
    bonds), plus the per-sublattice automaton MPO tensors themselves.
    Together these are literally the paper's "Hamiltonian terms outside the
    window, projected onto the reduced D-dimensional Hilbert space of the
    left/right block" (arXiv:1804.09163, Sec. V.1) -- an ordinary iDMRG
    growth environment already *is* that projection, nothing new needs to
    be solved for it. `window_boundary` is the 8-tuple
    `(HL, HL_bra, HL_ket, HL_mpo, HR, HR_bra, HR_ket, HR_mpo)`; `None`
    entries throughout mean no macro-iteration ever completed (e.g.
    maxiter's minimum of 2 was barely met). Consumed by
    `idmrg_window.py`'s heterogeneous-window builder, not meant to be used
    directly by ordinary callers."""

    def __init__(self, sites_uc, n_uc, U_list, e0, converged, niter_done,
                 state_overlap=None, local_superblock=None,
                 W_bulk=None, window_boundary=None):
        self.sites_uc = sites_uc
        self.n_uc = n_uc
        self.U_list = U_list
        self.e0 = e0
        self.converged = converged
        self.niter_done = niter_done
        self.state_overlap = state_overlap
        self.local_superblock = local_superblock
        self.W_bulk = W_bulk
        (self.env_HL, self.env_HL_bra, self.env_HL_ket, self.env_HL_mpo,
         self.env_HR, self.env_HR_bra, self.env_HR_ket, self.env_HR_mpo) = (
            window_boundary if window_boundary is not None else (None,) * 8)


def idmrg_ground_state(site_types, h_intra_op, h_inter_op, n_uc, maxm=30,
                        cutoff=1e-12, maxiter=200, etol=1e-10, niter=30,
                        verbose=False):
    """Run the growing/infinite-size DMRG algorithm to convergence (or
    `maxiter` macro-iterations) for a unit cell of `n_uc` sites (type codes
    `site_types`) and Hamiltonian given by `h_intra_op`/`h_inter_op` (plain
    MultiOperator.op-format term lists, 0-based site indices -- h_intra_op
    touching only sites 0..n_uc-1, h_inter_op touching at least one site
    n_uc..2*n_uc-1 -- see infinitechain.py's `_canonicalize_hamiltonian`
    for how a user-facing L/C/R-suffixed Hamiltonian is turned into this
    shape). Returns an IDMRGResult; `.e0` is the converged ground-state
    energy *per site* (not a total -- see this module's docstring)."""
    if n_uc > 2:
        # infinitechain.py's constructor already rejects n_uc>2 before
        # this function is ever reached through the public API -- this is
        # defense in depth for a direct caller of this module (bypassing
        # the wrapper): without it, n_uc=3 fails several micro-steps in
        # with an opaque "ValueError: axes don't match array" three stack
        # frames deep inside kernels.py's matvec, confirmed directly --
        # exactly the confusing-failure-mode this check exists to avoid.
        # See this module's own docstring / infinitechain.py's
        # __init__ for why n_uc>2 isn't supported yet.
        raise NotImplementedError(
            "idmrg_ground_state: n_uc>2 (got {}) is not supported yet -- "
            "see this module's docstring".format(n_uc))
    if maxiter < 2:
        # The energy *density* is a finite difference between two
        # consecutive macro-iterations' local ground-state eigenvalues
        # (see the loop below), so it is structurally undefined before a
        # second macro-iteration has run -- confirmed directly: maxiter=0
        # left IDMRGResult.e0=None with niter_done=1 (misreporting that an
        # iteration ran), and maxiter=1 likewise returned e0=None with no
        # error, only surfacing later as a confusing AttributeError deep
        # inside a subsequent vev()/correlator() call.
        raise ValueError(
            "idmrg_ground_state: maxiter must be >= 2 (energy density is "
            "a finite difference between two macro-iterations), got {}".format(maxiter))
    sites_uc, W_bulk = _build_automaton(h_intra_op, h_inter_op, site_types, n_uc)
    W_start0 = _project_channel(W_bulk[0], "left", 0)          # 0 = "S"
    W_endlast = _project_channel(W_bulk[n_uc - 1], "right", 1)  # 1 = "F"

    HL = HR = None
    HL_bra = HL_ket = None
    HR_bra = HR_ket = None
    # The Index currently serving as HL's/HR's own mpo axis -- tracked
    # explicitly rather than re-read off W_bulk[p] each time it's reused
    # (see the comment further below on why that reuse needs de-colliding).
    HL_mpo = HR_mpo = None
    energy = None
    prev_energy = None
    prev_density = None
    U_list = [None] * n_uc
    converged = False
    macro_iter = 0
    # prev_local[mstep]: the flattened local ground vector produced at this
    # unit-cell position by the *previous* macro-iteration, threaded back in
    # as the next macro-iteration's own Lanczos warm start -- see
    # _local_two_site_solve's own docstring for why. state_overlap tracks
    # how much the warm-started solve actually changed relative to what it
    # started from, refreshed every macro-iteration so the value returned
    # at the end always reflects the last macro-iteration actually run.
    prev_local = [None] * n_uc
    state_overlap = None
    # The very last micro-step's own local-solve ingredients (HL/HR
    # environments, MPO tensors, physical Indices, converged local ground
    # vector+energy) -- overwritten every micro-step, so it always holds the
    # *final* one actually run once the loop below exits (whether by
    # convergence or by exhausting maxiter). Kept so local_excitation_gap can
    # re-diagonalize this exact effective Hamiltonian for a second eigenpair
    # without rerunning the growing algorithm -- see IDMRGResult's own
    # docstring.
    last_superblock = None
    # Snapshot of (HL, HL_bra, HL_ket, HL_mpo, HR, HR_bra, HR_ket, HR_mpo) as
    # they stand *before* each macro-iteration's own mstep loop mutates them
    # -- overwritten every macro-iteration, so whatever it holds once the
    # loop exits (break on convergence, or maxiter exhaustion) is exactly
    # the environment entering the *last executed* macro-iteration. This is
    # deliberately not "the very latest HL/HR" (those are one mstep-loop
    # *ahead* -- they already absorbed the final macro-iteration's own
    # sites) -- it is instead consistent, by construction, with the
    # *returned* U_list's own edge bonds: U_list[0]'s own left bond is
    # exactly this snapshot's HL_ket (svd()'s left_inds always starts with
    # HL_ket, see _local_two_site_solve), and U_list[n_uc-1]'s own right
    # bond is exactly this snapshot's HR_ket, for the same reason on the
    # right. Used to cap a finite window with infinite boundary conditions
    # for real-time dynamical correlators (arXiv:1804.09163) -- see
    # IDMRGResult's own docstring.
    env_window_boundary = (None,) * 8

    for macro_iter in range(maxiter):
        env_window_boundary = (HL, HL_bra, HL_ket, HL_mpo, HR, HR_bra, HR_ket, HR_mpo)
        overlaps_this_iter = []
        for mstep in range(n_uc):
            p_L = mstep
            p_R = n_uc - 1 - mstep
            first_step = (macro_iter == 0 and mstep == 0)
            if first_step:
                W_pL, W_pR = W_start0, W_endlast
            else:
                # W_bulk[p] is the *same* tensor object every time
                # sublattice p comes up again (every macro-iteration past
                # the first), still carrying _build_automaton's own
                # boundary_idx identity on its left/right legs -- but
                # HL's/HR's *own* mpo axis was already re-minted fresh at
                # the end of the previous step (see below), so it must be
                # forced onto W_bulk[p]'s matching leg here by position,
                # not assumed to already agree. (Using the raw, un-forced
                # W_bulk[p] here instead measurably breaks the algorithm:
                # for n_uc=1, _build_automaton's single channel list makes
                # W_bulk[0]'s own left and right legs literally the same
                # Index object, which -- left un-broken -- makes HL's and
                # HR's own mpo axes collide onto that same object too,
                # corrupting every _local_two_site_solve from the second
                # macro-iteration on: confirmed directly by extracting the
                # dense effective Hamiltonian at macro-iteration 2 of a
                # uniform n_uc=1 Heisenberg chain and comparing its full
                # eigenvalue spectrum against exact 6-site ED -- completely
                # different spectra (not just truncation error), while
                # macro-iteration 1's spectrum matched ED to machine
                # precision, isolating the corruption to exactly this
                # reuse. The same collision risk exists for other small
                # n_uc, not just n_uc=1, since _build_automaton's
                # boundary_idx pool has only n_uc distinct objects total.)
                W_pL = _relabel_pos(W_bulk[p_L], 0, HL_mpo)
                W_pR = _relabel_pos(W_bulk[p_R], -1, HR_mpo)
            if p_L == p_R:
                W_pR = _fresh_physical_copy(W_pR)
            phys_L = _unprimed_site_index(W_pL)
            phys_R = _unprimed_site_index(W_pR)

            energy, U, S, V, evec0 = _local_two_site_solve(
                HL, HL_bra, HL_ket, W_pL, phys_L,
                W_pR, phys_R, HR, HR_bra, HR_ket,
                cutoff=cutoff, maxdim=maxm, niter=niter,
                x0_warm=prev_local[mstep])
            last_superblock = dict(
                HL=HL, HL_bra=HL_bra, HL_ket=HL_ket,
                W_pL=W_pL, phys_L=phys_L,
                W_pR=W_pR, phys_R=phys_R,
                HR=HR, HR_bra=HR_bra, HR_ket=HR_ket,
                evec0=evec0, energy=energy)
            if prev_local[mstep] is not None and prev_local[mstep].size == evec0.size:
                denom = np.linalg.norm(prev_local[mstep]) * np.linalg.norm(evec0)
                if denom > 0:
                    overlaps_this_iter.append(
                        abs(np.vdot(prev_local[mstep], evec0)) / denom)
            prev_local[mstep] = evec0

            left_ket_old, right_ket_old = HL_ket, HR_ket
            new_bond_u, new_bond_v = U.inds[-1], V.inds[0]

            # Fresh, independently-minted mpo-axis identities for the
            # *extend* step only (not used in the solve above, so this
            # doesn't disturb the W_pL/W_pR crossing bond the solve just
            # relied on) -- guarantees the new HL's and HR's own mpo axes
            # can never collide with each other, or with a future reuse of
            # W_bulk[p], regardless of how small n_uc is.
            new_HL_mpo = Index(W_pL.inds[-1].dim, tags="Link")
            W_pL_ext = _relabel_pos(W_pL, -1, new_HL_mpo)
            new_HR_mpo = Index(W_pR.inds[0].dim, tags="Link")
            W_pR_ext = _relabel_pos(W_pR, 0, new_HR_mpo)

            HL, HL_bra = _extend_HL(HL, HL_bra, W_pL_ext, U, left_ket_old, new_bond_u)
            HL_ket = new_bond_u
            HR, HR_bra = _extend_HR(HR, HR_bra, W_pR_ext, V, right_ket_old, new_bond_v)
            HR_ket = new_bond_v
            HL_mpo, HR_mpo = new_HL_mpo, new_HR_mpo

            U_list[p_L] = U

        state_overlap = min(overlaps_this_iter) if overlaps_this_iter else None
        density = ((energy - prev_energy) / (2 * n_uc)
                   if prev_energy is not None else None)
        if verbose:
            print("idmrg macro-iter {}: E={} density={} state_overlap={}".format(
                macro_iter, energy, density, state_overlap))
        if (density is not None and prev_density is not None
                and abs(density - prev_density) < etol):
            prev_energy, prev_density = energy, density
            converged = True
            break
        prev_energy, prev_density = energy, density

    if not converged and verbose:
        print("idmrg_ground_state: reached maxiter={} without converging "
              "to etol={} (last density change available)".format(maxiter, etol))

    return IDMRGResult(sites_uc, n_uc, U_list, prev_density, converged,
                        macro_iter + 1, state_overlap=state_overlap,
                        local_superblock=last_superblock,
                        W_bulk=W_bulk, window_boundary=env_window_boundary)


def local_excitation_gap(result, niter=200):
    """A cheap, cruder alternative to idmrg_excitations.excitation_energies/
    excitation_gap: the "local superblock gap" -- re-diagonalize the *same*
    2-site effective Hamiltonian the growing algorithm already solved for
    its ground state at the very last micro-step (result.local_superblock),
    but for its second-lowest eigenvalue instead, and return the difference.

    This is the direct infinite-chain analogue of a well-known finite-DMRG
    trick: at the last sweep, instead of only ever asking Lanczos for the
    lowest Ritz pair, ask the same local effective Hamiltonian for its two
    lowest and report their gap. Finite DMRG's own dedicated excited-state
    method (dmrg.py's overlap-penalty dmrg_excited/Chain.excited_states) is
    usually phrased as a Lagrange-multiplier/penalty problem -- minimize
    <psi|H|psi> subject to <psi|psi>=1 and <psi|psi_0>=0 -- because it
    re-sweeps the *whole chain* with psi_0 held fixed as an external MPS,
    and enforcing "orthogonal to a whole separate MPS at every local step"
    genuinely needs a penalty term threaded through the sweep. Here there is
    no separate sweep: the "psi_0" to stay orthogonal to is just the local
    ground vector already found, in the very same local Hilbert space, by
    the very same 2-site solve -- so the constraint can be enforced
    *exactly*, as a hard projector (deflation: P = I - |psi0><psi0|, the
    same orthogonal-complement idea idmrg_excitations.py's null-space gauge
    fixing V_L already uses, just against a single vector here instead of a
    whole tangent-space direction), rather than approximately via a finite
    penalty weight -- a Hermitian eigenproblem's constrained stationary
    points are exactly its *other* eigenvectors, so this is what a
    Lagrange-multiplier penalty converges to anyway as its weight -> inf,
    just obtained directly with nothing to tune.

    IMPORTANT CAVEAT -- this is a fundamentally cruder notion of "gap" than
    idmrg_excitations.excitation_energies/excitation_gap (the tangent-space
    quasiparticle ansatz): it has no momentum label at all (no e^{ikn}
    superposition over unit cells), so it reports a single number, not a
    dispersion, and it reuses HL/HR exactly as they converged for the
    *ground* state -- they are never allowed to relax/reoptimize for
    whatever the second eigenstate actually is. It is best understood as
    "the lowest-energy state the current, ground-state-optimized local
    ansatz has left over once the ground state itself is projected out",
    not a variationally-optimal excited state of the infinite chain. Use
    idmrg_excitations for a physically principled dispersion/gap when D=1
    applies; use this as a cheap, order-of-magnitude cross-check, or for
    the D>1 cases the tangent-space ansatz does not support yet.

    Returns a single float: the second-lowest local eigenvalue minus the
    ground-state eigenvalue (both real by Hermiticity of the local
    effective Hamiltonian)."""
    sb = result.local_superblock
    if sb is None:
        raise RuntimeError(
            "local_excitation_gap: this IDMRGResult has no stored local "
            "superblock (idmrg_ground_state must run at least one "
            "micro-step)")
    HL, HL_bra, HL_ket = sb["HL"], sb["HL_bra"], sb["HL_ket"]
    W_pL, phys_L = sb["W_pL"], sb["phys_L"]
    W_pR, phys_R = sb["W_pR"], sb["phys_R"]
    HR, HR_bra, HR_ket = sb["HR"], sb["HR_bra"], sb["HR_ket"]
    evec0, e0 = sb["evec0"], sb["energy"]

    order_in = ([HL_ket] if HL_ket is not None else []) + [phys_L, phys_R] + \
               ([HR_ket] if HR_ket is not None else [])
    shape = tuple(ind.dim for ind in order_in)
    s_L_out, s_R_out = phys_L.prime(1), phys_R.prime(1)
    order_out = ([HL_bra] if HL_bra is not None else []) + [s_L_out, s_R_out] + \
                ([HR_bra] if HR_bra is not None else [])
    pieces = [p for p in (HL, W_pL, W_pR, HR) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    dim = int(np.prod(shape))
    if dim < 2:
        raise RuntimeError(
            "local_excitation_gap: local Hilbert space has dimension {} -- "
            "too small to hold a state orthogonal to the ground "
            "state".format(dim))
    psi0 = evec0 / np.linalg.norm(evec0)

    def _deflate(v):
        return v - psi0 * np.vdot(psi0, v)

    if dim <= 3:
        # Same small-dim fallback as _local_two_site_solve: too small a
        # space for Lanczos to be meaningful, diagonalize directly.
        basis = np.eye(dim, dtype=complex)
        Hmat = np.column_stack([matvec(basis[:, k]) for k in range(dim)])
        w, _ = np.linalg.eigh((Hmat + Hmat.conj().T) / 2)
        e1 = w[1]
    else:
        def deflated_matvec(v):
            # P H P, P = I - |psi0><psi0| -- Hermitian since both P and the
            # underlying local Hamiltonian are (see this function's own
            # docstring), so _lanczos_ground_state's own Hermiticity
            # assumption still holds restricted to psi0's orthogonal
            # complement.
            return _deflate(matvec(_deflate(v)))
        rng = np.random.default_rng()
        v0 = _deflate(rng.standard_normal(dim) + 1j * rng.standard_normal(dim))
        e1_c, _ = _lanczos_ground_state(deflated_matvec, v0, niter=max(niter, 200))
        e1 = e1_c.real
    return float(e1 - e0)


def local_excitation_gap_windowed(result, h_intra_op, h_inter_op, site_types,
                                   n_uc, window=1, niter=300):
    """PROTOTYPE, not part of the stable public API yet -- exploring whether
    `local_excitation_gap`'s "frozen HL/HR" approximation is actually the
    dominant source of its error, by growing the local diagonalization block
    to `window` extra *free* physical sites on each side of the original 2
    (rather than adding Krylov vectors to the same frozen 2-site block, which
    `local_excitation_gap`'s own dim>3 branch already does via deflated
    Lanczos -- more Krylov vectors there cannot change the answer, since it's
    already an exact diagonalization of that specific effective Hamiltonian).
    `window=0` reduces to exactly the same effective Hamiltonian
    `local_excitation_gap` diagonalizes (used as an internal consistency
    check -- see below), just re-solving for its own ground state via Lanczos
    from scratch instead of reusing the growing algorithm's already-known
    `evec0` -- the two must agree to Lanczos precision, or something is wrong
    with the construction below.

    Only `n_uc=1` is supported: widening requires inserting extra copies of
    the periodic per-sublattice MPO tensor between HL/HR and the original 2
    sites, and for `n_uc=1` there is only one sublattice type, so every
    inserted copy is identical (`W_bulk[0]`) -- no bookkeeping needed for
    *which* sublattice position each extra site takes. `n_uc=2` would need
    that bookkeeping (alternating sublattice, and knowing which of the two
    the original phys_L/phys_R sit at) which result.local_superblock does not
    record and this prototype does not attempt.

    Construction: HL/HR (raw, unregularized, exactly as stored on
    result.local_superblock) are extended outward by `window` fresh copies of
    W_bulk[0] each (via `_fresh_physical_copy`, so every inserted site gets
    its own distinct physical Index -- colliding two sites onto the same
    Index is exactly the bug `_fresh_physical_copy`'s own docstring warns
    about) with freshly-minted Link Indices threading HL_mpo -> extra site 0
    -> extra site 1 -> ... -> (relabeled) W_pL's own left leg, mirroring
    exactly how `idmrg_ground_state`'s own loop threads W_pL/W_pR onto
    HL_mpo/HR_mpo via `_relabel_pos`. The extra sites' own physical legs are
    genuinely free (not contracted against any ground-state tensor) -- they
    become additional entries in the local eigenproblem's own `order_in`/
    `order_out`, exactly like phys_L/phys_R already are, so both the ground
    state AND its orthogonal-complement deflated first excited state are
    re-solved fresh within this larger space (unlike `local_excitation_gap`,
    which reuses the growing algorithm's already-known ground vector).

    IMPORTANT CAVEAT (measure, don't assume): HL/HR themselves are still
    exactly the growing algorithm's own converged, but *raw/unregularized*
    (not energy-density-subtracted) accumulator matrices -- widening the
    window does not relax them, it only lets the excitation spread across
    more *real* physical sites before hitting that frozen boundary. Whether
    this actually converges the gap toward the true answer as `window`
    grows, or by how much, is exactly the empirical question this prototype
    exists to answer -- not asserted here."""
    if n_uc != 1:
        raise NotImplementedError(
            "local_excitation_gap_windowed: only n_uc=1 is supported "
            "(got n_uc={}) -- see this function's own docstring".format(n_uc))
    sb = result.local_superblock
    if sb is None:
        raise RuntimeError(
            "local_excitation_gap_windowed: this IDMRGResult has no stored "
            "local superblock (idmrg_ground_state must run at least one "
            "micro-step)")
    HL, HL_bra, HL_ket = sb["HL"], sb["HL_bra"], sb["HL_ket"]
    W_pL, phys_L = sb["W_pL"], sb["phys_L"]
    W_pR, phys_R = sb["W_pR"], sb["phys_R"]
    HR, HR_bra, HR_ket = sb["HR"], sb["HR_bra"], sb["HR_ket"]
    if HL is None or HR is None:
        raise RuntimeError(
            "local_excitation_gap_windowed: cannot widen the window -- the "
            "stored superblock has an open (None) boundary on at least one "
            "side (run more iDMRG macro-iterations first)")

    def _mpo_axis(H, H_bra, H_ket):
        cands = [ind for ind in H.inds if ind != H_bra and ind != H_ket]
        if len(cands) != 1:
            raise RuntimeError(
                "local_excitation_gap_windowed: internal error -- expected "
                "exactly one non-bra/ket leg on the stored environment, "
                "found {}".format(len(cands)))
        return cands[0]

    HL_mpo = _mpo_axis(HL, HL_bra, HL_ket)
    HR_mpo = _mpo_axis(HR, HR_bra, HR_ket)

    _, W_bulk_fresh = _build_automaton(h_intra_op, h_inter_op, site_types, n_uc)
    W0 = W_bulk_fresh[0]

    extra_L = []
    bond = HL_mpo
    for _ in range(window):
        w = _relabel_pos(_fresh_physical_copy(W0), 0, bond)
        new_bond = Index(w.inds[-1].dim, tags="Link")
        w = _relabel_pos(w, -1, new_bond)
        extra_L.append(w)
        bond = new_bond
    W_pL_shifted = _relabel_pos(W_pL, 0, bond)

    extra_R = []
    bond = HR_mpo
    for _ in range(window):
        w = _relabel_pos(_fresh_physical_copy(W0), -1, bond)
        new_bond = Index(w.inds[0].dim, tags="Link")
        w = _relabel_pos(w, 0, new_bond)
        extra_R.append(w)
        bond = new_bond
    W_pR_shifted = _relabel_pos(W_pR, -1, bond)
    extra_R = list(reversed(extra_R))  # innermost (closest to phys_R) first

    phys_extra_L = [_unprimed_site_index(w) for w in extra_L]
    phys_extra_R = [_unprimed_site_index(w) for w in extra_R]
    bra_extra_L = [ind.prime(1) for ind in phys_extra_L]
    bra_extra_R = [ind.prime(1) for ind in phys_extra_R]

    order_in = ([HL_ket] if HL_ket is not None else []) + phys_extra_L + \
               [phys_L, phys_R] + phys_extra_R + \
               ([HR_ket] if HR_ket is not None else [])
    shape = tuple(ind.dim for ind in order_in)
    s_L_out, s_R_out = phys_L.prime(1), phys_R.prime(1)
    order_out = ([HL_bra] if HL_bra is not None else []) + bra_extra_L + \
                [s_L_out, s_R_out] + bra_extra_R + \
                ([HR_bra] if HR_bra is not None else [])
    pieces = [p for p in ([HL] + extra_L + [W_pL_shifted, W_pR_shifted] +
                           extra_R + [HR]) if p is not None]
    matvec = kernels.make_matvec(pieces, order_in, shape, order_out)

    dim = int(np.prod(shape))
    if dim < 2:
        raise RuntimeError(
            "local_excitation_gap_windowed: local Hilbert space has "
            "dimension {} -- too small to hold a state orthogonal to the "
            "ground state".format(dim))

    rng = np.random.default_rng()
    v0 = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    e0_c, psi0_c = _lanczos_ground_state(matvec, v0, niter=max(niter, 200))
    psi0 = psi0_c / np.linalg.norm(psi0_c)
    e0 = e0_c.real

    def _deflate(v):
        return v - psi0 * np.vdot(psi0, v)

    def deflated_matvec(v):
        return _deflate(matvec(_deflate(v)))

    v1 = _deflate(rng.standard_normal(dim) + 1j * rng.standard_normal(dim))
    e1_c, _ = _lanczos_ground_state(deflated_matvec, v1, niter=max(niter, 200))
    e1 = e1_c.real
    return float(e1 - e0)


# -- static correlators, via the standard infinite-MPS transfer-matrix
# formalism, operating on IDMRGResult.U_list (all plain NumPy from here on
# -- bond dimension/position is all that matters, not Index identity, so
# there is no need to carry ITensor Index bookkeeping into this part). ----

def _to_array_lpr(U):
    """U's array as (chi_left, d, chi_right) -- matches svd()'s own Index
    order (left_inds..., bond) exactly, so this is just U.array reshaped
    when U has no left bond at all (only possible for a genuinely
    first-ever site, which a *converged* unit cell's U_list never is)."""
    if len(U.inds) == 3:
        return U.array
    return U.array.reshape((1,) + U.array.shape)


def _transfer_matrices(U_list, n_uc, bra_list=None):
    """[E_p] for p=0..n_uc-1, each a (chi_l,chi_l,chi_r,chi_r) array (or
    (chi_l_ket,chi_l_bra,chi_r_ket,chi_r_bra) when `bra_list` differs from
    `U_list` -- see below): the doubled ket-times-conj(bra) transfer tensor
    for the tensor at sublattice p. Takes a plain tensor list rather than an
    IDMRGResult so apply_mpo's own canonicalization can reuse it on its own
    intermediate (not-yet-converged IDMRGResult) tensors too.

    `bra_list=None` (the default, used by every caller except
    `imps_overlap`) means the ordinary *self*-overlap transfer tensor,
    ket-times-conj(same tensor) -- what onsite_expectation/
    two_point_correlator/_canonicalize_periodic all need. Passing a
    *different* periodic tensor list computes the *mixed* transfer tensor
    between two distinct states instead -- what `imps_overlap` needs to
    compute <bra|ket> between two independently converged iMPS, which may
    not even share a bond dimension (hence the four-way (l_ket,l_bra,r_ket,
    r_bra) shape above, only square in the self-overlap case)."""
    if bra_list is None:
        bra_list = U_list
    Es = []
    for p in range(n_uc):
        A = _to_array_lpr(U_list[p])
        Bc = _to_array_lpr(bra_list[p])
        Es.append(np.einsum('lpr,LpR->lLrR', A, np.conj(Bc)))
    return Es


def _compose(E_a, E_b):
    """Transfer tensor for two sites in a row: E_a's own (r,R) legs
    contracted against E_b's (l,L) legs."""
    return np.einsum('lLrR,rRsS->lLsS', E_a, E_b)


def _apply_transfer(E4, rho):
    return np.einsum('lLrR,rR->lL', E4, rho)


# Relative-gap threshold below which two of the transfer matrix's leading
# eigenvalue magnitudes are treated as a genuine tie (see
# _dominant_right_fixed_point's own degeneracy check). The genuine tie
# this check exists to catch (two independently-converged, equally-
# normalized branches summed via imps_sum) lands at a ~1e-15 relative gap
# (exact to machine precision, confirmed directly on an up/down-polarized
# pair of product states -- idmrg.imps_sum's own docstring has the
# physical argument for why *any* two ordinary IDMRGResults tie here,
# always), so this threshold is set far tighter (1e-9) than that -- a
# huge margin, deliberately, rather than the loosest value that still
# passed known-legitimate cases: a *gapless/critical* single (non-summed)
# state's own top-two eigenvalue gap shrinks with growing maxm (finite-
# entanglement scaling -- confirmed directly on the uniform n_uc=1
# Heisenberg chain: gap 0.376 at maxm=40, 0.116 at maxm=60, a clear
# maxm-dependent trend, not noise), so a threshold merely "loose enough
# for the maxm values this repo's own tests happen to use today" would
# risk a *future* false positive on an ordinary, non-degenerate
# correlator call at a larger maxm a user might reasonably choose for
# better accuracy near criticality. 1e-9 keeps ~9 orders of magnitude of
# headroom below every legitimate gap measured so far while still
# rejecting the real tie (~1e-15) by 6 more orders of magnitude than
# needed -- a gapped (dimerized) n_uc=2 chain's top two eigenvalue
# magnitudes were (1.0, 0.093) and the gapless case above never dropped
# below a 0.116 gap up to maxm=60, both leaving comfortable margin above
# this threshold.
_DEGENERACY_RTOL = 1e-9


def _check_dominant_eigenvalue_nondegenerate(w, caller):
    """Raise RuntimeError if the two largest-magnitude entries of `w` are
    within `_DEGENERACY_RTOL` of each other, otherwise return `w`'s own
    descending-|.|  sort order -- shared by every "pick a single dominant
    eigenvalue/eigenvector" consumer in this module
    (`_dominant_right_fixed_point`, `_dominant_left_fixed_point`,
    `_dominant_eigenvalue_mixed`): a single dominant fixed point is not
    well-defined when the leading eigenvalue is (near-)degenerate --
    `np.argmax` would otherwise silently pick *one* arbitrary member of
    the tied eigenspace (confirmed directly: summing two independently
    converged, oppositely Sz-polarized `n_uc=1` product states via
    `idmrg.imps_sum` -- both individually normalized to eta=1, hence
    exactly tied here -- reliably reproduced only ONE of the two
    branches, with which one came out being sensitive to floating-point
    tie-breaking noise in `np.linalg.eig`, not a reproducible or even
    well-defined choice). This is not a rare edge case for
    `idmrg.imps_sum`'s own main use: any two separately-converged
    IDMRGResults are *always* individually normalized to eta=1 by SVD
    construction, so summing two genuinely different ground states hits
    exactly this tie essentially every time -- see `imps_sum`'s own
    docstring for the physical reason (a "cat state" superposition of two
    macroscopically distinct branches is not representable as a single
    injective/canonical periodic MPS) and why this is a documented,
    deliberate scope limit rather than a bug to route around here. `caller`
    is folded into the error message purely so it's attributable to
    whichever of this function's several callers actually raised, since
    `imps_sum` is only one of them (`apply_mpo`/`_canonicalize_periodic`,
    `onsite_expectation`, `two_point_correlator`, and `imps_overlap`'s own
    self-overlap terms all reach `_dominant_right_fixed_point` too, and a
    degenerate spectrum there need not have anything to do with
    `imps_sum` at all)."""
    order = np.argsort(-np.abs(w))
    if len(order) > 1 and abs(w[order[1]]) > (1 - _DEGENERACY_RTOL) * abs(w[order[0]]):
        raise RuntimeError(
            "idmrg ({}): the transfer matrix's dominant eigenvalue is "
            "(near-)degenerate (top two magnitudes {} and {}) -- a single "
            "dominant eigenvector is not well-defined here. This is most "
            "often seen when a state is, or derives from (e.g. via "
            "idmrg.imps_sum), a superposition of two macroscopically "
            "distinct branches with matched per-site norm -- see "
            "idmrg.imps_sum's own docstring for why that case is out of "
            "scope for this module's correlator machinery".format(
                caller, abs(w[order[0]]), abs(w[order[1]])))
    return order


def _dominant_right_fixed_point(Es):
    """The dominant right eigenvector of the full unit-cell transfer
    matrix T=E_0...E_{n_uc-1} (as a chi x chi density-matrix-like array,
    normalized to trace 1) and its eigenvalue -- should be close to 1 for
    a properly converged, correctly normalized (each U_list[p] isometric
    by SVD construction) infinite chain."""
    T4 = Es[0]
    for E in Es[1:]:
        T4 = _compose(T4, E)
    chi = T4.shape[0]
    if T4.shape != (chi, chi, chi, chi):
        # U_list[0]'s own left bond and U_list[n_uc-1]'s own right bond
        # are, in a truly periodic MPS, the *same* wraparound bond -- but
        # each is independently truncated by its own micro-step's SVD (up
        # to maxm, based on that cut's own entanglement spectrum), so the
        # growing algorithm's raw IDMRGResult.U_list offers no guarantee
        # they come out numerically equal, even at good convergence.
        # Confirmed directly: happens intermittently for some (maxm,
        # maxiter) combinations, not universally -- raise clearly here
        # instead of letting a generic reshape error surface deep inside.
        raise RuntimeError(
            "idmrg static correlators: the converged unit cell's wraparound "
            "bond dimension is inconsistent (U_list[0]'s left bond and "
            "U_list[-1]'s right bond differ, transfer tensor shape {}) -- "
            "try a different maxm/maxiter/etol combination for "
            "gs_energy()".format(T4.shape))
    Tmat = T4.reshape(chi * chi, chi * chi)
    w, v = np.linalg.eig(Tmat)
    order = _check_dominant_eigenvalue_nondegenerate(w, "_dominant_right_fixed_point")
    idx = order[0]
    eta = w[idx]
    rho = v[:, idx].reshape(chi, chi)
    rho = rho / np.trace(rho)
    return rho, eta


def _all_right_fixed_points(Es, n_uc):
    """rho_after[p] = the fixed-point "everything strictly after site p,
    wrapping back around" density matrix, for every sublattice position --
    obtained from one dominant-eigenvector computation (the p=n_uc-1 case)
    plus n_uc-1 cheap transfer-tensor applications, rather than a fresh
    eigenproblem per position."""
    rho_full, eta = _dominant_right_fixed_point(Es)
    rho_after = [None] * n_uc
    rho_after[n_uc - 1] = rho_full
    cur = rho_full
    for p in range(n_uc - 1, 0, -1):
        cur = _apply_transfer(Es[p], cur)
        cur = cur / np.trace(cur)
        rho_after[p - 1] = cur
    return rho_after, eta


def _op_transfer(sites_uc, U_list, p, opname):
    """Transfer tensor for sublattice p with operator `opname` inserted
    (applied to the ket side only, as <bra|opname|ket> per site)."""
    M = sites_uc.site_type(p + 1).matrix(opname)
    A = _to_array_lpr(U_list[p])
    Aop = np.einsum('io,lir->lor', M, A)
    return np.einsum('lpr,LpR->lLrR', Aop, np.conj(A))


def onsite_expectation(result, opname, p):
    """<opname> at sublattice p (0..n_uc-1) of the converged infinite
    chain."""
    Es = _transfer_matrices(result.U_list, result.n_uc)
    rho_after, _eta = _all_right_fixed_points(Es, result.n_uc)
    E_op = _op_transfer(result.sites_uc, result.U_list, p, opname)
    val = _apply_transfer(E_op, rho_after[p])
    return complex(np.trace(val))


def two_point_correlator(result, opname_i, p_i, opname_j, r):
    """<opname_i(site p_i) opname_j(site p_i + r)> of the converged
    infinite chain, r measured in physical sites (r=0: both operators at
    the same site, composed as opname_i then opname_j, i.e. matrix product
    M_i @ M_j in the (in,out) convention of sites/base.py -- matching how
    `SxC[i]*SxC[i]`-style same-site products read left to right elsewhere
    in dmrgpy). r must be >= 0; use r and swap the operator order for the
    mirror separation if r<0 is needed."""
    if r < 0:
        raise ValueError("two_point_correlator: r must be >= 0")
    n_uc = result.n_uc
    Es = _transfer_matrices(result.U_list, n_uc)
    rho_after, _eta = _all_right_fixed_points(Es, n_uc)

    if r == 0:
        M = (result.sites_uc.site_type(p_i + 1).matrix(opname_i)
             @ result.sites_uc.site_type(p_i + 1).matrix(opname_j))
        A = _to_array_lpr(result.U_list[p_i])
        Aop = np.einsum('io,lir->lor', M, A)
        E4 = np.einsum('lpr,LpR->lLrR', Aop, np.conj(A))
        val = _apply_transfer(E4, rho_after[p_i])
        return complex(np.trace(val))

    p_j = (p_i + r) % n_uc
    running = _op_transfer(result.sites_uc, result.U_list, p_i, opname_i)
    for k in range(1, r):
        p = (p_i + k) % n_uc
        running = _compose(running, Es[p])
    running = _compose(running, _op_transfer(
        result.sites_uc, result.U_list, p_j, opname_j))
    val = _apply_transfer(running, rho_after[p_j])
    return complex(np.trace(val))


def _dominant_eigenvalue_mixed(Es):
    """Dominant eigenvalue of the full unit-cell *mixed* transfer matrix
    T=E_0...E_{n_uc-1}, generalizing `_dominant_right_fixed_point`'s own
    eigen-solve to the case where the ket and bra tensor lists don't share
    a single common bond dimension -- exactly the case `imps_overlap` needs
    (the two states being overlapped may have been converged/truncated to
    different maxm). Unlike `_dominant_right_fixed_point`, this does not
    require T4's four legs to all share one size, only that the *ket*'s own
    wraparound bond (legs 0 and 2) and the *bra*'s own wraparound bond (legs
    1 and 3) are each individually self-consistent -- the two need not
    match each other. Returns only the eigenvalue: `imps_overlap` has no
    use for the fixed-point vector itself, unlike the correlator machinery
    that consumes `_dominant_right_fixed_point`'s own rho. Also runs the
    same degeneracy check as `_dominant_right_fixed_point`/
    `_dominant_left_fixed_point` (see
    `_check_dominant_eigenvalue_nondegenerate`'s own docstring) -- a tied
    top eigenvalue here means "which branch dominates" is not well-defined
    for the cross term either, the same underlying issue, just on the
    mixed rather than self-overlap transfer matrix. Confirmed this does
    not spuriously trip on realistic `imps_overlap` inputs: two
    independently-converged ordinary ground states' own mixed transfer
    spectra (both a genuinely orthogonal case, where the mixed matrix
    reduces to a single trivial eigenvalue with nothing to tie against,
    and a genuinely overlapping case) leave a wide margin, and
    `imps_overlap(result, apply_mpo(result, W_identity))`'s own
    gauge-comparison case (a non-trivial bond dimension, unlike the
    orthogonal case above) does too."""
    T4 = Es[0]
    for E in Es[1:]:
        T4 = _compose(T4, E)
    chi_ket, chi_bra = T4.shape[0], T4.shape[1]
    if T4.shape != (chi_ket, chi_bra, chi_ket, chi_bra):
        raise RuntimeError(
            "idmrg imps_overlap: a state's own wraparound bond dimension "
            "is inconsistent (mixed transfer tensor shape {}) -- same "
            "failure mode _dominant_right_fixed_point already guards "
            "against for a single state, see its own comment".format(T4.shape))
    Tmat = T4.reshape(chi_ket * chi_bra, chi_ket * chi_bra)
    w = np.linalg.eigvals(Tmat)
    order = _check_dominant_eigenvalue_nondegenerate(w, "_dominant_eigenvalue_mixed")
    return w[order[0]]


def imps_overlap(result_a, result_b, normalize=True):
    """Per-unit-cell overlap <result_b|result_a> between two converged
    infinite MPS (`result_a`/`result_b`: IDMRGResult and/or PeriodicMPS, any
    combination -- both are accepted anywhere that only reads
    .sites_uc/.n_uc/.U_list, matching onsite_expectation/
    two_point_correlator's own duck typing).

    A literal <psi_b|psi_a> over an infinite chain is not, in general, a
    finite number: it factorizes into one power of the dominant *mixed*
    transfer-matrix eigenvalue eta_ab (see `_dominant_eigenvalue_mixed`)
    per unit cell, so <psi_b|psi_a> ~ eta_ab^N is 0 or infinite as N (the
    number of unit cells) -> infinity, unless |eta_ab|==1 exactly. The
    physically meaningful, always-finite quantity is instead the *per-site
    fidelity*

        eta_ab / sqrt(eta_aa * eta_bb)

    where eta_aa/eta_bb are `result_a`/`result_b`'s own self-overlap
    eigenvalues (exactly 1 for a properly left-canonical U_list -- which is
    what both idmrg_ground_state and apply_mpo always return -- but
    computed here rather than assumed, so this also works on a hand-built
    or un-normalized PeriodicMPS). This per-site fidelity has magnitude 1
    iff `result_a` and `result_b` represent the same physical state (any
    gauge, any normalization convention on the raw tensors), and magnitude
    < 1 otherwise -- a transfer-matrix Cauchy-Schwarz bound that holds
    because every U_list tensor produced by this module is isometric/
    left-canonical. This is the standard iMPS notion of state overlap
    (e.g. used for fidelity susceptibility across a phase transition, or to
    check that `apply_mpo` with an identity-equivalent operator reproduced
    the original state up to gauge -- see test_imps_overlap_* in
    tests/test_infinite_chain.py). Pass normalize=False for the raw,
    un-normalized eta_ab instead -- mainly a diagnostic, analogous to
    apply_mpo's own returned `eta` (which is exactly this function's
    would-be eta_aa, applied to apply_mpo's output against itself).

    Requires both states to share the same n_uc and, at every sublattice
    position, the same local physical dimension (same Hilbert space) --
    raises ValueError otherwise. Bond dimensions need not match between the
    two states (this is exactly why the mixed-transfer eigen-solve below
    uses `_dominant_eigenvalue_mixed`, not `_dominant_right_fixed_point`,
    which assumes a single self-consistent chi)."""
    n_uc = result_a.n_uc
    if result_b.n_uc != n_uc:
        raise ValueError(
            "imps_overlap: unit-cell size mismatch (result_a.n_uc={}, "
            "result_b.n_uc={})".format(n_uc, result_b.n_uc))
    dims_a = [_to_array_lpr(result_a.U_list[p]).shape[1] for p in range(n_uc)]
    dims_b = [_to_array_lpr(result_b.U_list[p]).shape[1] for p in range(n_uc)]
    if dims_a != dims_b:
        raise ValueError(
            "imps_overlap: physical dimension mismatch per sublattice "
            "(result_a={}, result_b={})".format(dims_a, dims_b))

    Es_ab = _transfer_matrices(result_a.U_list, n_uc, bra_list=result_b.U_list)
    eta_ab = _dominant_eigenvalue_mixed(Es_ab)
    if not normalize:
        return complex(eta_ab)

    _, eta_aa = _dominant_right_fixed_point(_transfer_matrices(result_a.U_list, n_uc))
    _, eta_bb = _dominant_right_fixed_point(_transfer_matrices(result_b.U_list, n_uc))
    return complex(eta_ab) / np.sqrt(complex(eta_aa) * complex(eta_bb))


# == Applying a (bounded) MPO to the converged iMPS ================
#
# apply_mpo(result, W_bulk) is the infinite-chain analogue of
# mpsalgebra.applyMPO(): contract a periodic MPO onto the converged unit
# cell (grow_by_mpo), then re-canonicalize/truncate the grown bond
# dimension back down (_canonicalize_periodic) via the standard two-sided
# fixed-point procedure (Orus & Vidal, "Infinite time-evolving block
# decimation algorithm beyond unitary evolution", PRB 78, 155117 (2008)) --
# the same construction iTEBD uses after every non-unitary gate
# application, generalized here from a single-site unit cell to n_uc sites.
#
# SCOPE: W_bulk must represent a *bounded* (non-extensive) periodic
# operator -- the same tensor reused at every unit cell, with no
# unconditional "keep accumulating forever" self-loop. `_build_periodic_mpo`
# (the Hamiltonian's own automaton builder, used by idmrg_ground_state
# above) is deliberately the *other* kind: its "F" channel has an
# unconditional identity self-loop specifically so a *finite,
# boundary-terminated* growing chain can sum an unbounded number of
# repeats (see _build_periodic_mpo's own docstring). Feeding that
# Hamiltonian automaton directly into apply_mpo below -- a boundary-less
# periodic contraction -- does NOT compute "H|psi>": traced around a ring
# with no boundary projection, the accumulator channel *multiplies* itself
# every unit cell instead of *summing* (worked through by hand on a
# trivial onsite-only automaton with no bond terms: going around N times
# literally computes (Id+onsite)^N, not sum_i onsite_i). apply_mpo is
# therefore scoped to genuinely bounded/local periodic operators --
# single-site products, gates tiled once per unit cell, symmetry
# operators, and the like -- built directly (see the docstrings below) or
# via `_build_periodic_mpo` fed a Hamiltonian-shaped term list that is
# genuinely finite-range *and* has no term wrapping past one adjacent
# cell's own automaton reach; using the Hamiltonian's own full automaton
# is out of scope for this function.


def _primed_site_index(W):
    """The primed ("out") Site-tagged Index of an MPO tensor W -- mirrors
    idmrg_ground_state's own `_unprimed_site_index`, used by grow_by_mpo to
    label the still-primed physical leg of a freshly MPO-grown tensor
    before apply_mpo noPrime()s it back to a plain ket leg."""
    return next(ind for ind in W.inds if ind.hastags("Site") and ind.plev == 1)


def _local_grow(W_p, A_p):
    """W_p applied to A_p at one unit-cell site: contract the shared
    physical leg and Kronecker-merge (left(W),left(A)) and
    (right(W),right(A)) into the two output bond axes -- done via plain
    positional NumPy array manipulation (W_p's own axis order is always
    (left,phys-in,phys-out,right), A_p's is always (left,phys,right), by
    construction, see _build_periodic_mpo/_to_array_lpr) rather than
    ITensor.__mul__'s Index-identity-based matching.

    Index-identity matching would silently misbehave here whenever W_p's
    own left and right legs happen to be the *same* Index object --
    guaranteed for n_uc=1 (see _build_periodic_mpo's own docstring on its
    boundary_idx pool having only n_uc distinct Index objects, reused
    verbatim for a repeating channel type) -- exactly the same collision
    idmrg_ground_state's own `_project_channel`/`_relabel_pos` already
    have to sidestep positionally for identical reasons.

    Returns a plain (chi_W_left*chi_A_left, d, chi_W_right*chi_A_right)
    NumPy array -- the physical leg is W_p's own "out" (still primed)
    convention, grow_by_mpo noPrime()s it once the full periodic chain is
    assembled."""
    w_arr = W_p.array
    a_arr = _to_array_lpr(A_p)
    Lw, d, _d2, Rw = w_arr.shape
    La, _d3, Ra = a_arr.shape
    prod = np.einsum('lsor,msn->lmorn', w_arr, a_arr)
    return prod.reshape((Lw * La, d, Rw * Ra))


def grow_by_mpo(W_bulk, U_list, n_uc):
    """Contract MPO tensor W_bulk[p] against ket tensor U_list[p] at every
    unit-cell site (via `_local_grow`), Kronecker-merging each of the
    n_uc *cuts* (including the wraparound one, between site n_uc-1 and
    site 0) into one freshly-minted combined Link Index -- the same
    per-site "zip together, don't compress yet" construction
    mpsalgebra.py's `_apply_chain` uses for a finite chain
    (mpsalgebra.py:206-251), just over a periodic index range instead of a
    chain with two open ends.

    W_bulk[p]'s physical (unprimed) Index must be the *same* object as
    U_list[p]'s own physical Index, so `_local_grow`'s contraction lands on
    exactly that leg -- true whenever W_bulk was built by
    `_build_periodic_mpo` against the *same* sites_uc U_list's own physical
    legs came from (see apply_mpo's own docstring).

    Returns B[0..n_uc-1]: rank-3 ITensors (left Link, phys -- still primed,
    W's own "out" convention, see apply_mpo which noPrime()s it -- right
    Link), *not yet* canonicalized/truncated (see `_canonicalize_periodic`)."""
    if len(W_bulk) != n_uc or len(U_list) != n_uc:
        raise ValueError(
            "grow_by_mpo: expected {} unit-cell sites, got W_bulk={}, "
            "U_list={}".format(n_uc, len(W_bulk), len(U_list)))
    raw = [_local_grow(W_bulk[p], U_list[p]) for p in range(n_uc)]
    combined_links = [Index(raw[p].shape[0], tags="Link") for p in range(n_uc)]
    for p in range(n_uc):
        right_dim = raw[p].shape[-1]
        expected = combined_links[(p + 1) % n_uc].dim
        if right_dim != expected:
            raise RuntimeError(
                "grow_by_mpo: cut dimension mismatch at site {} (right "
                "dim {}, next site's left dim {}) -- W_bulk and U_list "
                "must both be periodic with period n_uc={}".format(
                    p, right_dim, expected, n_uc))
    B = []
    for p in range(n_uc):
        phys_out = _primed_site_index(W_bulk[p])
        left, right = combined_links[p], combined_links[(p + 1) % n_uc]
        B.append(ITensor((left, phys_out, right), raw[p]))
    return B


def _apply_transfer_from_left(E4, rho):
    """Mirror of `_apply_transfer`: contracts rho against E4's *left*
    (l,L) legs instead of its right (r,R) ones -- used to propagate a
    dominant *left* fixed point forward through a site, the mirror image
    of how `_apply_transfer` propagates a right fixed point backward."""
    return np.einsum('lL,lLrR->rR', rho, E4)


def _dominant_left_fixed_point(Es):
    """The dominant LEFT eigenvector of the full unit-cell transfer matrix
    T=E_0...E_{n_uc-1} -- i.e. a vector rho_L with rho_L . T = eta * rho_L
    (equivalently, the right eigenvector of T's own transpose) -- at the
    *same* cut `_dominant_right_fixed_point` itself resolves (the
    wraparound bond, "before site 0"/"after site n_uc-1"). Mirrors
    `_dominant_right_fixed_point`'s own construction and
    shape-consistency check exactly, with one added transpose."""
    T4 = Es[0]
    for E in Es[1:]:
        T4 = _compose(T4, E)
    chi = T4.shape[0]
    if T4.shape != (chi, chi, chi, chi):
        raise RuntimeError(
            "idmrg apply_mpo: the periodic unit cell's wraparound bond "
            "dimension is inconsistent (transfer tensor shape {}) -- same "
            "failure mode _dominant_right_fixed_point already guards "
            "against, see its own comment".format(T4.shape))
    Tmat = T4.reshape(chi * chi, chi * chi)
    w, v = np.linalg.eig(Tmat.T)
    # Degeneracy check included here too (not just relying on
    # _canonicalize_periodic always calling _dominant_right_fixed_point
    # first, see _check_dominant_eigenvalue_nondegenerate's own
    # docstring) -- Tmat.T has the same eigenvalues as Tmat, so this is
    # the same tie condition, just independently guarded here in case a
    # future caller reaches this function without going through the
    # right-fixed-point check first.
    order = _check_dominant_eigenvalue_nondegenerate(w, "_dominant_left_fixed_point")
    idx = order[0]
    eta = w[idx]
    rho = v[:, idx].reshape(chi, chi)
    rho = rho / np.trace(rho)
    return rho, eta


def _all_left_fixed_points(Es, n_uc):
    """rho_before[p] = the fixed-point "everything strictly before site p,
    wrapping back around" density matrix, for every sublattice position --
    mirrors `_all_right_fixed_points`, propagating *forward* (rho_before[0]
    is the dominant left fixed point itself; rho_before[p] for p>0 comes
    from pushing it through sites 0..p-1) instead of backward.

    Also returns `scales`: `_all_right_fixed_points`' own per-step
    `cur = cur / np.trace(cur)` renormalization (needed there, and copied
    here, purely to avoid numerical blow-up/decay across many propagation
    steps) is harmless for `_all_right_fixed_points`' own callers -- the
    correlator code only ever uses one rho_after[p] at a time, and
    `_canonicalize_periodic`'s G_right[p] construction is provably
    invariant to rho_R_before[p]'s own overall scale (the S_p factor it
    divides by scales the same way, see `_canonicalize_periodic`'s own
    comment). G_left[p], by contrast, is *linear* in rho_L_before[p]'s own
    scale (no compensating S factor) -- so silently renormalizing away the
    *relative* scale between different cuts' rho_before here would corrupt
    `_canonicalize_periodic`'s left-canonicality by exactly that dropped
    per-cut ratio. Confirmed directly: for n_uc=2, this was the actual bug
    behind a 0.875/(1/0.875) Gram-matrix deviation at the two sites (their
    product is exactly 1, i.e. a real, non-unitary rescale, not numerical
    noise) that survived an unrelated (and, it turned out, inert -- see
    that function's own comment) rho_R normalization fix; n_uc=1 has no
    propagation step at all (scales=[1]) which is why it was never
    affected. `scales[p]` is the accumulated product of exactly the trace
    factors divided out through step p -- multiplying `rho_before[p]` by
    `scales[p]` undoes the renormalization exactly (by linearity of the
    transfer map, verified algebraically: if cur_normalized[p] =
    cur_natural[p]/scales[p], propagating one more (linear) step and
    renormalizing again preserves that relationship with
    scales[p+1]=scales[p]*this step's own divided-out trace)."""
    rho_full, eta = _dominant_left_fixed_point(Es)
    rho_before = [None] * n_uc
    rho_before[0] = rho_full
    scales = [1.0] * n_uc
    cur = rho_full
    scale = 1.0
    for p in range(0, n_uc - 1):
        cur = _apply_transfer_from_left(Es[p], cur)
        step_trace = np.trace(cur)
        cur = cur / step_trace
        scale = scale * step_trace
        scales[p + 1] = scale
        rho_before[p + 1] = cur
    return rho_before, eta, scales


def _psd_sqrt_factor(rho, rel_floor=1e-12):
    """X (k x N, k <= N) such that rho ~= X^H X (Hermitized square root
    via eigh, dropping eigenvalues at or below `rel_floor` times the
    largest one -- this covers both clipping small-negative numerical
    noise, since 0 <= rel_floor*max, and dropping small-but-positive ones)
    -- factors the dominant left/right transfer-matrix fixed points ahead
    of `_canonicalize_periodic`'s SVD-based truncation/regauging step. Not
    `svd.py`'s `eigh_truncate`: that also truncates to a caller-chosen
    cutoff/maxdim (the *final* truncation `_canonicalize_periodic` must
    only apply once, jointly, from the SVD of X @ Y, not here) and only
    returns eigenvectors, not a full sqrt factor.

    Dropping near-zero eigenvalues here (not just clipping negative ones)
    is required, not cosmetic: a raw environment fixed point can be
    enormously ill-conditioned (confirmed directly on a genuinely
    bond-growing apply_mpo case -- eigenvalues spanning ~1e-11 to ~1,
    condition number ~1e10) whenever the *raw* grown bond dimension (an
    artifact of grow_by_mpo's Kronecker product, before any truncation)
    vastly exceeds the state's actual entanglement at that cut. Keeping
    those near-singular directions in a *square*, full-rank X/Y produces
    an M_p = X_p @ Y_p whose own SVD is numerically garbage in exactly
    those directions -- confirmed directly: without this floor, a genuine
    (non-uniform-bond-dimension) grow-then-truncate case left one site's
    Gram matrix off from Identity by O(100), not a small residual, while
    the same run's *other* site (whose own environment was well-
    conditioned) was fine -- narrowing the cause to exactly this."""
    herm = (rho + rho.conj().T) / 2
    evals, evecs = np.linalg.eigh(herm)
    floor = rel_floor * evals[-1] if evals.size and evals[-1] > 0 else 0.0
    keep = evals > max(floor, 0.0)
    return np.sqrt(evals[keep])[:, None] * evecs[:, keep].conj().T


def _canonicalize_periodic(B_list, n_uc, cutoff, maxdim):
    """The standard two-sided fixed-point infinite-MPS canonicalization/
    compression procedure (see this module's "Applying a (bounded) MPO"
    section docstring above for the reference). B_list: raw (generally
    non-canonical) periodic tensors, already noPrime()d, one per
    unit-cell site, n_uc cuts total (including the wraparound).

    At each cut p, factor the dominant left/right fixed points
    rho_L_before[p] = X_p^H X_p, rho_R_before[p] = Y_p Y_p^H
    (`_psd_sqrt_factor`), SVD M_p = X_p @ Y_p = U_p S_p V_p^H, truncate
    (`svd.py`'s `_truncate`), and build the *asymmetric* gauge
    G_left[p] = U_p^H X_p (no S factor at all), G_right[p] =
    Y_p V_p^H S_p^-1 (the *full* inverse, not a square root) -- inserting
    G_right[p] @ G_left[p] at cut p still resolves (up to the truncation
    just performed) the identity there (S cancels: U^H X Y V^H S^-1 =
    U^H (U S V^H) V^H S^-1 = S S^-1 = Identity_keep on the kept subspace,
    same argument as a symmetric S^-1/2 split, just redistributed), but
    -- unlike a symmetric S^-1/2 split, which produces Vidal Gamma tensors
    still needing a separate bond-weight Lambda=S layer threaded between
    them -- putting the *entire* S^-1 on one side directly produces the
    plain left-canonical form this module's IDMRGResult.U_list already
    uses (each tensor alone isometric, no separate bond-weight layer),
    with no further Lambda-absorption step needed. (Confirmed numerically
    on the identity-MPO round-trip case before committing to this
    formula: the symmetric-S^-1/2-plus-separate-Lambda-absorption
    construction a first derivation produced looked plausible by analogy
    to Vidal's canonical form but was measurably wrong -- gram matrix
    deviation from Identity ~1, not a small residual; this asymmetric
    form reproduces Identity to ~1e-14.)

    Returns (U_list_new, eta): U_list_new is left-canonical (verified
    internally, raising if it isn't within tolerance -- a real bug would
    otherwise silently corrupt every downstream correlator computation on
    the result); eta is the new state's own self-overlap transfer
    eigenvalue (a norm diagnostic -- apply_mpo does not renormalize, so
    this is not necessarily close to 1)."""
    Es = _transfer_matrices(B_list, n_uc)
    rho_after, eta_R = _all_right_fixed_points(Es, n_uc)
    rho_R_before = [rho_after[(p - 1) % n_uc] for p in range(n_uc)]
    rho_L_before, eta_L, scales = _all_left_fixed_points(Es, n_uc)

    if abs(eta_L - eta_R) > 1e-6 * max(1.0, abs(eta_R)):
        raise RuntimeError(
            "idmrg apply_mpo: dominant left/right transfer eigenvalues "
            "disagree (eta_L={}, eta_R={}) -- same transfer matrix, only "
            "the eigenvector side differs, so a mismatch signals a "
            "genuine inconsistency (e.g. a near-degenerate transfer "
            "spectrum) rather than benign numerical noise".format(eta_L, eta_R))

    # Undo _all_left_fixed_points' own per-cut renormalization -- see that
    # function's own comment for why G_left[p] (unlike G_right[p]) is not
    # invariant to it, and why this is the fix (not, as an earlier,
    # confirmed-inert attempt assumed, a mismatch between rho_L and rho_R).
    #
    # Also transpose: _apply_transfer_from_left(E4, rho)[r,R] pairs its
    # *first* input index with E4's 'l' leg (the *unconjugated*/ket copy)
    # and its second with 'L' (the conjugated/bra copy) -- see
    # _transfer_matrices' own E4 = einsum('lpr,LpR->lLrR', A, conj(A)).
    # That makes _all_left_fixed_points' own rho_before a (ket,bra)-ordered
    # matrix, the *transpose* of the usual (bra,ket) density-matrix
    # convention the isometry identity below needs (the same convention
    # rho_R_before/rho_after already use, consumed correctly as-is).
    # Confirmed directly, the hard way: an isometry-target identity derived
    # by hand (G_right[p]^H @ rho_L_before[p] @ G_right[p] == Identity_keep,
    # with rho_L_before[p] appearing via
    # sum_{b,b'} rho_L_before[0][b,b'] conj(B[b,..]) B[b',..] == rho_L_before[1])
    # checked out in the identity/algebra but the *numeric* Gram matrix of
    # the resulting tensor was still off by O(10-1000) (not the earlier
    # bugs' signatures, and unaffected by cutoff/maxdim/psd-factor-rank
    # choices) until this transpose was added -- isolated by computing
    # that same bracketed sum directly and finding it equal to
    # conj(rho_L_before[1]) (== rho_L_before[1].T for a Hermitian matrix),
    # not rho_L_before[1] itself.
    rho_L_before = [(rho_L_before[p] * scales[p]).T for p in range(n_uc)]

    G_left, G_right = [None] * n_uc, [None] * n_uc
    for p in range(n_uc):
        X_p = _psd_sqrt_factor(rho_L_before[p])
        Y_p = _psd_sqrt_factor(rho_R_before[p]).conj().T
        M_p = X_p @ Y_p
        U_p, S_p, Vh_p = np.linalg.svd(M_p, full_matrices=False)
        keep, _discarded = _truncate(S_p, cutoff, maxdim, mindim=1)
        U_p, S_p, Vh_p = U_p[:, :keep], S_p[:keep], Vh_p[:keep, :]
        G_left[p] = U_p.conj().T @ X_p
        G_right[p] = (Y_p @ Vh_p.conj().T) * (1.0 / S_p)[None, :]

    U_list_new = []
    for p in range(n_uc):
        arr = np.einsum('ab,bpc,cd->apd', G_left[p], B_list[p].array,
                         G_right[(p + 1) % n_uc])
        left = Index(arr.shape[0], tags="Link")
        right = Index(arr.shape[2], tags="Link")
        phys = next(ind for ind in B_list[p].inds if ind.hastags("Site"))
        U_list_new.append(ITensor((left, phys, right), arr))

    for p in range(n_uc):
        arr = U_list_new[p].array
        gram = np.einsum('apc,apd->cd', arr.conj(), arr)
        ident = np.eye(gram.shape[0])
        # 1e-4, not machine precision: X_p/Y_p's own conditioning degrades
        # with how far the *raw* grown bond dimension (chi_A*chi_W, before
        # any truncation) exceeds the state's real entanglement at that
        # cut -- confirmed directly on a deliberately extreme stress case
        # (an already near-maxdim-saturated 16-dim bond grown by a
        # bond-4 gate with maxdim=None, so almost nothing is discarded):
        # a genuine, essentially unavoidable ~1e-10 condition number left
        # a ~1e-5 to 1e-3 residual depending on cutoff/maxdim, even though
        # the exact same construction reproduces Identity to ~1e-13 on
        # every less pathological case tried (n_uc in (1,2), identity and
        # unitary chi_W=1 MPOs, and this same gate at a smaller starting
        # bond dimension). 1e-4 is loose enough to absorb that, while
        # still being orders of magnitude below every *actual* bug this
        # check caught during development (deviations of 1, 50-1000s, or
        # more -- see this function's own docstring history).
        if not np.allclose(gram, ident, atol=1e-4):
            raise RuntimeError(
                "idmrg apply_mpo: canonicalization produced a non-left-"
                "canonical tensor at sublattice {} (max deviation from "
                "Identity: {}) -- indicates a bug in the gauge-fixing "
                "construction, not a benign numerical issue".format(
                    p, np.max(np.abs(gram - ident))))

    return U_list_new, eta_R


class PeriodicMPS:
    """A periodic iMPS with no ground-state-specific bookkeeping -- same
    shape as IDMRGResult (sites_uc, n_uc, U_list) minus e0/converged/
    niter_done, which have no meaning for apply_mpo's output.
    onsite_expectation/two_point_correlator only ever read
    .sites_uc/.n_uc/.U_list off their `result` argument, so they accept a
    PeriodicMPS directly, no changes needed there. `eta` is apply_mpo's own
    diagnostic (see `_canonicalize_periodic`'s docstring)."""

    def __init__(self, sites_uc, n_uc, U_list, eta):
        self.sites_uc = sites_uc
        self.n_uc = n_uc
        self.U_list = U_list
        self.eta = eta


def apply_mpo(result, W_bulk, cutoff=1e-12, maxdim=None):
    """Apply a periodic MPO to the converged iMPS `result` (an IDMRGResult
    or PeriodicMPS), returning a new PeriodicMPS representing W|psi> up to
    (cutoff, maxdim) truncation -- the infinite-chain analogue of
    `mpsalgebra.applyMPO`. See this module's "Applying a (bounded) MPO to
    the converged iMPS" section docstring above for the *scope*
    restriction on W_bulk (bounded/local periodic operators only -- not
    the Hamiltonian's own unbounded automaton) and the algorithm
    (`grow_by_mpo` + `_canonicalize_periodic`)."""
    B = grow_by_mpo(W_bulk, result.U_list, result.n_uc)
    B = [_t_noPrime(b, "Site") for b in B]
    U_list_new, eta = _canonicalize_periodic(B, result.n_uc, cutoff, maxdim)
    return PeriodicMPS(result.sites_uc, result.n_uc, U_list_new, eta)


# == Summing two converged iMPS ====================================
#
# imps_sum(result_a, result_b) is the periodic-chain analogue of
# mpsalgebra.sum(): the standard, exact MPS-addition construction
# (block-diagonal in the bond space, `_periodic_direct_sum`), followed by
# the same `_canonicalize_periodic` re-gauging/truncation `apply_mpo`
# already uses. Unlike mpsalgebra.sum's finite chain -- whose two open
# ends are dimension-1, so the block-diagonal bulk construction plus a
# plain concatenation at the boundaries reproduces a literal Hilbert-space
# vector sum exactly -- a periodic chain has no open end: every cut is a
# "bulk" cut, so the direct sum is block-diagonal *everywhere*, including
# the wraparound bond.
#
# PHYSICAL SCOPE -- read before treating this as a drop-in infinite
# analogue of mpsalgebra.sum: tiled to the thermodynamic limit, a
# block-diagonal unit cell does not represent |psi_a>+|psi_b> the way the
# finite construction represents a literal vector sum. Tracing it around N
# unit cells gives eta_a^N + eta_b^N (eta_a/eta_b: result_a/result_b's own
# self-overlap transfer eigenvalues, always ~1 exactly for any
# IDMRGResult, by left-canonical SVD construction -- see
# `_local_two_site_solve`'s own comment), so:
#
# - Two branches with a genuine norm mismatch (|eta_a| != |eta_b|, only
#   possible for a hand-built or already-rescaled PeriodicMPS -- not two
#   ordinary IDMRGResults, which are always both ~1) have a well-posed,
#   non-degenerate combined transfer matrix: the smaller-norm branch is
#   exponentially suppressed as N -> infinity, so `_canonicalize_periodic`
#   correctly (not a bug -- the true infinite-volume answer) returns a
#   state numerically indistinguishable from the larger-norm branch alone.
# - Two branches with matched |eta| -- which is the *common* case, since
#   summing two separately-converged ground states (e.g. two
#   symmetry-related solutions of the same translationally-invariant
#   Hamiltonian, the most natural reason to want this operation at all)
#   always ties exactly -- produce a genuinely degenerate dominant
#   eigenspace. Representing that "cat state" (a superposition of two
#   macroscopically distinct branches) as a single injective/canonical
#   periodic MPS is not possible with the fixed-point machinery this
#   module's correlator functions rely on everywhere else, so
#   `_dominant_right_fixed_point`'s own degeneracy check (see its
#   docstring, calibrated directly against real single-state spectra)
#   reliably raises RuntimeError in this case instead of silently
#   returning an arbitrary, wrong single branch -- confirmed directly:
#   before that check existed, summing two independently-converged,
#   oppositely Sz-polarized product states (both eta=1, orthogonal)
#   silently collapsed to just one of the two branches, chosen by
#   floating-point tie-breaking noise inside np.linalg.eig, not even
#   reproducibly. Correctly representing this common, physically
#   meaningful case (e.g. via the standard thermodynamic-limit local-
#   observable identity <O> -> (<O>_a+<O>_b)/2 for two exactly orthogonal,
#   equally-weighted branches) needs correlator machinery this module
#   does not have yet -- a documented, deliberate scope limit, in the
#   same spirit as apply_mpo's own "bounded operators only" restriction
#   above, not something silently gotten wrong.


def _periodic_direct_sum(U_list_a, U_list_b, n_uc):
    """Raw (uncanonicalized) block-diagonal direct sum of two periodic
    tensor lists, one cut per unit-cell site (n_uc cuts total, including
    the wraparound) -- the periodic-chain analogue of mpsalgebra.sum's own
    per-site block-diagonal placement (mpsalgebra.py:143-203), applied at
    *every* cut instead of leaving the two end cuts as a plain
    concatenation (there is no open end here to concatenate along instead
    -- see this module's "Summing two converged iMPS" section docstring
    above for the physical meaning of this once tiled to the
    thermodynamic limit).

    U_list_a[p]/U_list_b[p] must share the same physical dimension at
    every sublattice position p (checked explicitly below, redundantly
    with `imps_sum`'s own pre-check, since this function is also called
    directly in tests/other internal code -- a mismatched physical
    dimension would otherwise silently broadcast the smaller physical
    leg across the larger one instead of raising, confirmed directly by
    hand-checking the equivalent NumPy assignment pattern). The returned
    tensor's own physical Index at position p is U_list_a[p]'s own --
    U_list_b[p]'s physical Index need not be the same object, only the
    same dimension, since this function only ever matches legs by
    position (`_to_array_lpr`), never by Index identity.

    Also requires U_list_a and U_list_b to each independently be
    consistent periodic chains (site p's own right bond dimension must
    equal site (p+1)%n_uc's own left bond dimension, within each list
    separately) -- mirrors `grow_by_mpo`'s own analogous cut-dimension
    check (this function's closest sibling in this module) rather than
    letting a mismatch surface as an opaque NumPy broadcast error deep in
    the array-assignment loop below."""
    arrs_a = [_to_array_lpr(U_list_a[p]) for p in range(n_uc)]
    arrs_b = [_to_array_lpr(U_list_b[p]) for p in range(n_uc)]
    for p in range(n_uc):
        if arrs_a[p].shape[1] != arrs_b[p].shape[1]:
            raise ValueError(
                "_periodic_direct_sum: physical dimension mismatch at "
                "sublattice {} (U_list_a={}, U_list_b={})".format(
                    p, arrs_a[p].shape[1], arrs_b[p].shape[1]))
        next_p = (p + 1) % n_uc
        if arrs_a[p].shape[2] != arrs_a[next_p].shape[0]:
            raise RuntimeError(
                "_periodic_direct_sum: U_list_a's own bond dimension is "
                "inconsistent at cut {} (site {}'s right dim {}, site "
                "{}'s left dim {}) -- U_list_a must itself be a "
                "consistent periodic chain".format(
                    p, p, arrs_a[p].shape[2], next_p, arrs_a[next_p].shape[0]))
        if arrs_b[p].shape[2] != arrs_b[next_p].shape[0]:
            raise RuntimeError(
                "_periodic_direct_sum: U_list_b's own bond dimension is "
                "inconsistent at cut {} (site {}'s right dim {}, site "
                "{}'s left dim {}) -- U_list_b must itself be a "
                "consistent periodic chain".format(
                    p, p, arrs_b[p].shape[2], next_p, arrs_b[next_p].shape[0]))
    combined_links = [Index(arrs_a[p].shape[0] + arrs_b[p].shape[0], tags="Link")
                       for p in range(n_uc)]
    out = []
    for p in range(n_uc):
        Aa, Ab = arrs_a[p], arrs_b[p]
        La, d, Ra = Aa.shape
        left, right = combined_links[p], combined_links[(p + 1) % n_uc]
        arr = np.zeros((left.dim, d, right.dim), dtype=complex)
        arr[:La, :, :Ra] = Aa
        arr[La:, :, Ra:] = Ab
        phys = next(ind for ind in U_list_a[p].inds if ind.hastags("Site"))
        out.append(ITensor((left, phys, right), arr))
    return out


def imps_sum(result_a, result_b, cutoff=1e-12, maxdim=None):
    """Direct sum of two converged infinite MPS (`result_a`/`result_b`:
    IDMRGResult and/or PeriodicMPS, any combination -- same duck typing as
    `imps_overlap`/`apply_mpo`), returning a new PeriodicMPS:
    `_periodic_direct_sum`'s raw block-diagonal construction,
    re-canonicalized/truncated via the same `_canonicalize_periodic`
    two-sided fixed-point procedure `apply_mpo` already uses. See this
    module's "Summing two converged iMPS" section docstring above for what
    this "+" actually means once tiled to the thermodynamic limit -- in
    particular, why summing two *ordinary* IDMRGResults (always
    individually normalized to eta=1) reliably raises RuntimeError rather
    than silently returning one arbitrary branch.

    Requires both states to share the same n_uc and, at every sublattice
    position, the same local physical dimension (mirrors imps_overlap's
    own check) -- raises ValueError otherwise. Bond dimensions need not
    match."""
    n_uc = result_a.n_uc
    if result_b.n_uc != n_uc:
        raise ValueError(
            "imps_sum: unit-cell size mismatch (result_a.n_uc={}, "
            "result_b.n_uc={})".format(n_uc, result_b.n_uc))
    dims_a = [_to_array_lpr(result_a.U_list[p]).shape[1] for p in range(n_uc)]
    dims_b = [_to_array_lpr(result_b.U_list[p]).shape[1] for p in range(n_uc)]
    if dims_a != dims_b:
        raise ValueError(
            "imps_sum: physical dimension mismatch per sublattice "
            "(result_a={}, result_b={})".format(dims_a, dims_b))
    raw = _periodic_direct_sum(result_a.U_list, result_b.U_list, n_uc)
    U_list_new, eta = _canonicalize_periodic(raw, n_uc, cutoff, maxdim)
    return PeriodicMPS(result_a.sites_uc, n_uc, U_list_new, eta)
