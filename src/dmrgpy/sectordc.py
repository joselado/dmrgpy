"""Sector-resolved dynamical correlators: `submode="SECTOR"`.

The Lehmann representation of a zero-temperature dynamical correlator,

    C_AB(w) = sum_n <0|A|n> <n|B|0> * L_delta(w - (E_n - E_0)),

summed not over excited states of the whole Hilbert space (which is what
`submode="EX"`/dcex.py does) but over the lowest eigenstates of the *one
quantum-number sector that B|0> actually lands in*.

That sector is not a choice: `B` shifts every conserved charge by a fixed
amount, so `<n|B|0>` vanishes identically unless `|n>` carries the ground
state's charges plus `charge(B)`. For an electron spectral function on a
chain with N particles, `B = Cdag_i` means the intermediate states are the
N+1 eigenstates and nothing else; the hole part lives entirely in N-1. For
a spin chain with an Sz=0 ground state, `S+` reaches Sz=+1, `S-` reaches
Sz=-1 (Sz=+/-2 in the 2*Sz units the sector API uses throughout), and `Sz`
stays put.

Why bother, next to KPM/CVM/TD:

  * Every intermediate state is a separately converged DMRG eigenstate of
    a *smaller* Hilbert space. Symmetry, not an overlap penalty, is what
    keeps them out of the ground state and out of each other's charge
    channels -- they live in a different sector, so the search cannot
    collapse back into |0>. Within the target sector they are still found
    by the same overlap-penalty method `submode="EX"` uses, only in a much
    smaller space, and the usual fluctuation warning still applies to
    them. There is no broadening/resolution tradeoff, no
    expansion order, no time window and no Fourier artifact -- the poles
    come out at their converged energies with their exact weights, and
    `return_poles=True` hands them back as poles rather than as a curve.
  * The answer is resolved by quantum-number channel, which no other
    submode here is.

What it costs: only the lowest `nex` states of the sector are found, so
only the low-energy end of the response is described. That truncation is
the method's single approximation and it is measured rather than assumed
-- see the `captured` diagnostic below, which reports what fraction of the
spectral weight of B|0> the returned poles account for. This is the right
tool for a few sharp low-lying excitations and the wrong one for a broad
continuum; that is exactly the opposite tradeoff to KPM.

How it is possible at all: `promote_mps` rebases a wavefunction onto the
chain's *original* dense site indices, which `set_conserved_sector` never
rebuilds (it replaces `sites_`, not `dense_sites_`). So states solved in
different sectors of one live session can be promoted onto common indices
and contracted against each other -- which is what a cross-sector matrix
element like <n^{N+1}|Cdag_i|0^N> is. That is regression-tested
independently of this module, in
tests/test_sector_promotion.py::test_states_promoted_from_different_sectors_are_comparable.

Global phases need no fixing. Each converged |n> carries an arbitrary,
genuinely run-dependent phase (v3 starts from a random in-sector state),
but it cancels exactly in the product <0|A|n><n|B|0> -- the bra
contributes exp(-i theta_n), the ket exp(+i theta_n) -- as long as both
factors are taken from the same promoted handle, which they are. Only
inside a degenerate multiplet are the *individual* pole weights
basis-dependent; the multiplet sum is not.

itensor_version=3 and itensor_version="python" only: no other backend has
quantum numbers at all, and there is deliberately no ED fallback, since
falling back would answer with the global excited states -- a different
calculation wearing the same name.
"""

import warnings

import numpy as np

from . import operatornames
from .multioperatortk import charge as chargetk

# How many (reference sector, target sector, nex, ...) solves to keep
# alive at once. Three covers the spin structure factor's three channels
# and the spectral function's two, which are the sweeps that matter;
# beyond that the oldest entry goes, since each pins a clone chain.
_STATE_CACHE_SIZE = 4


def _resolve_operators(self, name, i, j):
    """(A,B) as MultiOperators. Normally already normalized by
    Many_Body_Chain.get_dynamical_correlator, but normalized again here so
    the module also works when called directly (dcex.py does the same)."""
    if name is None:
        raise ValueError(
            "get_dynamical_correlator(submode=\"SECTOR\"): name (the two-"
            "operator tuple, e.g. name=(fc.C[0],fc.Cdag[0])) must be given")
    name = operatornames.str2MO(self, name, i=i, j=j,
            require_symbolic_for="submode='SECTOR'")
    return name[0], name[1]


def _check_backend(self):
    """Refuse, rather than silently compute something else, on every
    backend/Hamiltonian this method does not apply to.

    This has to be explicit here. mode.py's own sector guard keys on
    `self.conserved_sector`, and the sector this method uses lives on an
    internal clone -- the caller's chain typically has none set, so that
    guard never fires for a submode="SECTOR" call.
    """
    if self.itensor_version not in (3, "python"):
        raise NotImplementedError(
            "get_dynamical_correlator(submode=\"SECTOR\") needs conserved-"
            "sector support, which exists for itensor_version=3 and "
            "itensor_version=\"python\" only (this chain uses %s). Switch "
            "with setup_cpp(version=3) or setup_python(), or use "
            "submode=\"KPM\"/\"TD\"/\"CVM\", which work on every backend."
            % repr(self.itensor_version))
    if self.itensor_version == 3 and self.ns < 3:
        raise NotImplementedError(
            "get_dynamical_correlator(submode=\"SECTOR\") on a %d-site "
            "chain: ITensor v3's two-site DMRG aborts below 3 sites, and "
            "the ED fallback cannot target a sector. Use "
            "itensor_version=\"python\" for a chain this short."
            % self.ns)
    if not self.is_hermitian(self.hamiltonian):
        raise NotImplementedError(
            "get_dynamical_correlator(submode=\"SECTOR\") assumes a "
            "Hermitian Hamiltonian: it sums a Lehmann representation over "
            "orthogonal eigenstates of one sector, which a non-Hermitian "
            "Hamiltonian does not provide. Use submode=\"KPM\"/\"CVM\", "
            "which have non-Hermitian implementations.")


def _quantum_numbers(self, clone, sector, conserve):
    """Which quantities to conserve, and the reference sector's targets.

    Three ways the reference sector gets fixed, in order of directness:
    the caller passes it, the chain already has one set (the natural
    workflow -- solve at fixed filling, then ask for the spectral
    function there), or it is measured on an unconstrained ground state
    and rounded. The measurement is a fallback rather than the default
    path because a genuinely charge-superposed ground state has no sector
    interpretation at all, and rounding one would invent an answer; that
    case warns loudly instead.
    """
    offered = chargetk.offered_qn_names(self)
    if not offered:
        raise ValueError(
            "get_dynamical_correlator(submode=\"SECTOR\"): no site of this "
            "chain carries a conserved quantum number, so there are no "
            "sectors to resolve (parafermion chains are the case). Use "
            "another submode.")
    if sector is not None:
        sector = {str(k): int(v) for k, v in dict(sector).items()}
    chain_sector = dict(self.conserved_sector or {})

    if conserve is not None:
        names = tuple(conserve)
    elif sector:
        names = tuple(sector.keys())
    elif chain_sector:
        names = tuple(chain_sector.keys())
    else:
        names = chargetk.conserved_qn_names(self, self.hamiltonian, offered)
        if not names:
            raise ValueError(
                "get_dynamical_correlator(submode=\"SECTOR\"): this "
                "Hamiltonian conserves none of the quantum numbers its "
                "sites offer (%s), so no sector is well defined. Pass "
                "conserve=/sector= explicitly if a quantity is conserved "
                "in a way this check cannot see."
                % ", ".join(offered))
    for name in names:
        if name not in offered:
            raise ValueError(
                "get_dynamical_correlator(submode=\"SECTOR\"): no site of "
                "this chain carries the quantum number %s (this chain "
                "offers %s)" % (repr(name), ", ".join(offered) or "none"))

    targets = {}
    for name in names:
        if sector and name in sector:
            targets[name] = sector[name]
        elif name in chain_sector:
            targets[name] = chain_sector[name]
        else:
            targets[name] = _measure_charge(self, clone, name)
    return tuple(names), targets


def _measure_charge(self, clone, name):
    """The reference sector's target for one quantity, measured on the
    unconstrained ground state of the clone (never on the caller's own
    chain, whose cached state this must not disturb)."""
    op = chargetk.charge_operator(self, name)
    if op is None:
        raise ValueError(
            "get_dynamical_correlator(submode=\"SECTOR\"): cannot work out "
            "which %s sector to start from -- this chain exposes no total-"
            "%s operator to measure it with. Say which sector explicitly, "
            "with set_conserved_sector(%s=...) on the chain or sector={%s: "
            "...} on this call." % (name, name, name, repr(name)))
    q = complex(clone.vev(op)).real
    target = int(np.round(q))
    if abs(q - target) > 1e-3:
        warnings.warn(
            "get_dynamical_correlator(submode=\"SECTOR\"): the "
            "unconstrained ground state has <%s> = %.6g, which is not an "
            "integer -- it is a superposition of sectors, and a "
            "sector-resolved correlator of it is not well defined. "
            "Proceeding from the nearest sector, %s=%d; pass sector= "
            "explicitly if that is not what you meant."
            % (name, q, name, target))
    return target


def _cap_nex(self, nex, target, quiet=False):
    """`nex` clipped to the number of states the target sector actually
    holds. Without this the default (20) overshoots small sectors -- the
    N+/-1 sectors of a short chain, an Sz=+2 sector of a six-site spin
    chain -- and the overlap-penalty search then spends most of the
    runtime hunting states that do not exist, warns about every one of
    them, and has its junk quietly discarded downstream. Slow, alarming,
    and it still returns a plausible-looking spectrum."""
    dim = chargetk.sector_dimension(self, target)
    if dim is not None and dim <= 0:
        raise ValueError(
            "get_dynamical_correlator(submode=\"SECTOR\"): the target "
            "sector %s is empty -- no state of this chain carries those "
            "quantum numbers, so the correlator is identically zero."
            % _sector_str(target))
    if dim is not None and nex > dim:
        if not quiet:
            warnings.warn(
                "get_dynamical_correlator(submode=\"SECTOR\"): nex=%d was "
                "reduced to %d, the dimension of the target sector %s -- "
                "that is every state it holds, so the Lehmann sum is "
                "complete rather than truncated."
                % (nex, dim, _sector_str(target)))
        nex = dim
    return max(int(nex), 1)


def _sector_str(sector):
    return ", ".join("%s=%d" % kv for kv in sorted(sector.items())) or "(none)"


def _sector_str_any(sector):
    """Same, for the one-or-many target sectors a decomposed correlator
    has (sector_lehmann turns target_sector into a list of them)."""
    if isinstance(sector, list):
        return " and ".join(_sector_str(s) for s in sector)
    return _sector_str(sector)



def _sector_states(self, clone, reference, target, nex, **kwargs):
    """The promoted ground state of `reference` and the `nex` lowest
    promoted states of `target`, all living on one session's dense site
    indices, cached on the chain.

    The cache is what makes a site sweep affordable. Nothing in the
    expensive half of this calculation depends on *which* operators are
    being correlated: one set of sector solves serves the entire matrix of
    (i,j) pairs, and a spectral function over every site of a chain would
    otherwise redo the same two DMRG solves once per site. It is keyed on
    everything that changes the states -- the two sectors, nex, the
    Hamiltonian's identity, the DMRG accuracy settings and the backend --
    and dropped by restart()/set_hamiltonian(), same as dcex.py's
    excited-state cache.

    The cached clone is kept alive with the states, deliberately: an MPS
    handle only means anything inside the session that minted it (Index
    identity is per session), so caching wavefunctions without their chain
    would cache something unusable.
    """
    key = (tuple(sorted(reference.items())), tuple(sorted(target.items())),
           int(nex), id(self.hamiltonian), self.itensor_version,
           self.maxm, self.nsweeps, self.cutoff, self.noise,
           tuple(sorted(kwargs.items(), key=lambda kv: kv[0])))
    # A dict, not a single slot: a spectral function alternates between
    # two target sectors (N+1 and N-1) as it sweeps sites, and a one-entry
    # cache would be evicted by every single call -- measured, that turned
    # a 6-site sweep into 6 full pairs of sector solves.
    cache = getattr(self, "_sector_states_cache", None)
    if not isinstance(cache, dict):
        cache = {}
        self._sector_states_cache = cache
    if key in cache:
        return cache[key]

    # -- the reference sector: ground state, promoted to dense indices
    clone.set_conserved_sector(**reference)
    e0 = clone.gs_energy()
    wf0 = clone.promote_mps(clone.get_gs())

    # -- the target sector: the nex lowest states, promoted likewise
    clone.set_conserved_sector(**target)
    energies, wfs = clone.get_excited_states(n=nex, purify=False, **kwargs)
    wfs = [clone.promote_mps(w) for w in wfs]

    # -- out of sector mode: only now can a charge-changing operator be
    # built at all (sector_terms rejects one inside a sector), and every
    # handle promoted above still matches the session's dense sites.
    clone.set_conserved_sector()

    out = (clone, wf0, wfs, [float(np.real(e)) for e in energies],
           float(np.real(e0)))
    if len(cache) >= _STATE_CACHE_SIZE:  # bounded: each entry pins a whole
        cache.pop(next(iter(cache)))     # clone chain (and its session)
    cache[key] = out
    return out


def _prepare(self, name, i, j, sector, conserve):
    """Everything the pipeline needs before any sector is touched: the two
    operators, the clone every sector switch happens on, and the reference
    sector."""
    A, B = _resolve_operators(self, name, i, j)
    _check_backend(self)
    clone = self.copy()  # every sector switch happens here, never on the
    # caller's chain: the pipeline ends outside the sector, and that must
    # not be a side effect the caller sees.
    names, reference = _quantum_numbers(self, clone, sector, conserve)
    return A, B, clone, names, reference


def _indefinite_charge_error(which, names):
    return ValueError(
        "get_dynamical_correlator(submode=\"SECTOR\"): the %s operator "
        "carries no definite %s charge, and cannot be split into parts "
        "that do, so B|gs> is not confined to any set of sectors and this "
        "method does not apply. Write the correlator in operators of "
        "definite charge (ladder operators for spin, Cdag/C for "
        "fermions), or use submode=\"KPM\"/\"TD\", which have no such "
        "restriction." % (which, "/".join(names)))


def _homogeneous_pairs(self, A, B, names):
    """The correlator as a sum of contributions with definite charge:
    `[(A_part, B_part, charge_of_B_part), ...]`.

    One pair when both operators already carry a definite charge, which is
    the common case and behaves exactly as if this function did not exist.
    Otherwise each operator is split into its charge components
    (`charge_components`) and the components are paired up: only
    charge(A_part) = -charge(B_part) survives, since every other product
    has <0|A_part|n><n|B_part|0> = 0 term by term. That is what makes
    name=(Sx,Sx) work -- Sx = (S+ + S-)/2, so it becomes the sum of an
    Sz+2 and an Sz-2 contribution, which is also the physically right
    decomposition.
    """
    dqA = chargetk.multioperator_charge(self, A, names)
    dqB = chargetk.multioperator_charge(self, B, names)
    if dqA is not None and dqB is not None:
        if any(a + b != 0 for a, b in zip(dqA, dqB)):
            raise ValueError(
                "get_dynamical_correlator(submode=\"SECTOR\"): the two "
                "operators change the conserved quantum numbers by %s and "
                "%s, which do not cancel -- <0|A|n><n|B|0> is then zero "
                "for every state and the correlator vanishes identically "
                "by symmetry. Check the operator order (the first operator "
                "is applied to the bra, the second to the ket)."
                % (_charge_str(names, dqA), _charge_str(names, dqB)))
        return [(A, B, dqB)]
    parts_A = chargetk.charge_components(self, A, names)
    parts_B = chargetk.charge_components(self, B, names)
    if parts_A is None:
        raise _indefinite_charge_error("first", names)
    if parts_B is None:
        raise _indefinite_charge_error("second", names)
    pairs = []
    for q in sorted(parts_B):
        opposite = tuple(-x for x in q)
        if opposite in parts_A:
            pairs.append((parts_A[opposite], parts_B[q], q))
    if not pairs:
        raise ValueError(
            "get_dynamical_correlator(submode=\"SECTOR\"): no charge "
            "component of the first operator (%s) cancels one of the "
            "second (%s), so <0|A|n><n|B|0> vanishes for every state and "
            "the correlator is identically zero by symmetry."
            % (", ".join(_charge_str(names, q) for q in sorted(parts_A)),
               ", ".join(_charge_str(names, q) for q in sorted(parts_B))))
    return pairs


def _poles_for_pair(self, clone, names, reference, A, B, dqB, nex,
                    quiet=False, **kwargs):
    """One charge-homogeneous contribution: the poles of <0|A|n><n|B|0>
    over the states of the single sector B|gs> lands in."""
    target = {n: reference[n] + d for n, d in zip(names, dqB)}
    nex = _cap_nex(self, nex, target, quiet=quiet)
    clone, wf0, wfs, energies, e0 = _sector_states(
        self, clone, reference, target, nex, **kwargs)

    norm0 = complex(clone.overlap(wf0, wf0)).real
    poles, weights, mel = [], [], []
    for e, wf in zip(energies, wfs):
        nn = complex(clone.overlap(wf, wf)).real
        scale = 1.0 / np.sqrt(nn * norm0)
        a = complex(clone.aMb(wf0, A, wf)) * scale
        b = complex(clone.aMb(wf, B, wf0)) * scale
        poles.append(float(np.real(e)) - float(np.real(e0)))
        weights.append(a * b)
        mel.append(b)

    total = complex(clone.aMb(wf0, B.get_dagger() * B, wf0)).real / norm0
    captured = float(sum(abs(x) ** 2 for x in mel) / total) if total > 1e-12 else 1.0
    info = {"reference_sector": reference, "target_sector": target,
            "quantum_numbers": names,
            "charge_A": tuple(-x for x in dqB), "charge_B": dqB,
            "nex": len(poles), "e0": float(np.real(e0)),
            "energies_absolute": [float(np.real(e)) for e in energies],
            "matrix_elements": np.array(mel), "captured": captured,
            "sum_rule_total": total}
    return np.array(poles), np.array(weights), info


def sector_poles(self, name=None, i=0, j=0, nex=20, sector=None,
                 conserve=None, quiet=False, **kwargs):
    """The Lehmann poles of a sector-resolved dynamical correlator over
    ONE sector: `(energies, weights, info)`.

    `energies` are excitation energies E_n - E_0 (they can be negative --
    adding a particle below the chemical potential lowers the energy);
    `weights` are the complex products <0|A|n><n|B|0>; `info` carries the
    sectors used, the absolute energies, the per-state matrix elements
    <n|B|0> and the captured spectral weight.

    Both operators must carry a definite charge here, so that a single
    target sector exists and `info` can describe it. `sector_lehmann`
    below is the general entry point, which splits a charge-indefinite
    operator (Sx, Sy) into parts first and sums their contributions;
    `dynamical_correlator` uses that one. This function stays the
    single-sector primitive, since a caller that wants the matrix elements
    themselves -- to build S(q,w) out of them, say -- needs them attached
    to one identified sector.
    """
    A, B, clone, names, reference = _prepare(self, name, i, j, sector,
                                             conserve)
    pairs = _homogeneous_pairs(self, A, B, names)
    if len(pairs) > 1:
        raise ValueError(
            "sector_poles: this operator pair spans %d charge channels "
            "(%s) and so has no single target sector. Use "
            "sector_lehmann()/get_dynamical_correlator(submode=\"SECTOR\"), "
            "which sums them, or pass operators of definite charge."
            % (len(pairs), ", ".join(_charge_str(names, q)
                                     for _, _, q in pairs)))
    Aq, Bq, dqB = pairs[0]
    out = _poles_for_pair(self, clone, names, reference, Aq, Bq, dqB, nex,
                          quiet=quiet, **kwargs)
    self._sector_dc_info = out[2]  # readable afterwards without changing
    # the (x,y) return shape every other submode has
    return out


def sector_lehmann(self, name=None, i=0, j=0, nex=20, sector=None,
                   conserve=None, quiet=False, **kwargs):
    """The full sector-resolved Lehmann poles: `(energies, weights, info)`
    with the contributions of every charge channel concatenated.

    Identical to `sector_poles` for the usual case of two operators of
    definite charge. When the operators are charge-indefinite -- name=
    (Sx,Sx) being the case that matters, since that is how essentially
    every spin correlator in dmrgpy's examples is written -- they are
    split into charge components and the channels whose charges cancel are
    summed, so the answer is the same correlator, obtained from the Sz+2
    and Sz-2 sectors instead of from one.
    """
    A, B, clone, names, reference = _prepare(self, name, i, j, sector,
                                             conserve)
    pairs = _homogeneous_pairs(self, A, B, names)
    poles, weights, channels = [], [], []
    for Aq, Bq, dqB in pairs:
        e, w, ch = _poles_for_pair(self, clone, names, reference, Aq, Bq,
                                   dqB, nex, quiet=quiet, **kwargs)
        poles.append(e)
        weights.append(w)
        channels.append(ch)
    info = dict(channels[0])
    if len(channels) > 1:
        # A single "the" target sector no longer means anything, and the
        # captured weight is the weight-averaged one across channels.
        total = sum(c["sum_rule_total"] for c in channels)
        info["target_sector"] = [c["target_sector"] for c in channels]
        info["captured"] = (sum(c["captured"] * c["sum_rule_total"]
                                for c in channels) / total
                            if total > 1e-12 else 1.0)
        info["sum_rule_total"] = total
        info["nex"] = sum(c["nex"] for c in channels)
        info["charge_A"] = [c["charge_A"] for c in channels]
        info["charge_B"] = [c["charge_B"] for c in channels]
        info["matrix_elements"] = np.concatenate(
            [c["matrix_elements"] for c in channels])
    info["channels"] = channels
    self._sector_dc_info = info
    return np.concatenate(poles), np.concatenate(weights), info


def _charge_str(names, charge):
    parts = ["%s=%+d" % (n, c) for n, c in zip(names, charge) if c != 0]
    return ", ".join(parts) or "(no charge)"


def dynamical_correlator(self, name=None, i=0, j=0, nex=20, delta=2e-2,
                         es=np.linspace(-1.0, 10.0, 1000), sector=None,
                         conserve=None, return_poles=False,
                         weight_tol=0.9, quiet=False, **kwargs):
    """Dynamical correlator from the eigenstates of one quantum-number
    sector -- see this module's docstring.

    Returns `(es, correlator)`, in exactly the convention every other
    submode uses (dcex.py's `1j*(advanced-retarded)/(2 pi)`, i.e. the
    complex array whose real part is the spectral function the ED
    reference `mode="ED", submode="ED"` returns). With
    `return_poles=True`, returns `(es, correlator, info)`, `info` also
    holding the raw poles -- which is the output this method exists for,
    since its poles are converged eigenvalues rather than features of a
    broadened curve.
    """
    poles, weights, info = sector_lehmann(
        self, name=name, i=i, j=j, nex=nex, sector=sector,
        conserve=conserve, quiet=quiet, **kwargs)
    if not quiet and info["captured"] > 1.0 + 1e-3:
        # More weight than B|gs> has: the states are not as orthogonal as
        # the sum assumes, so some of it is being counted twice. The
        # overlap-penalty search is still what separates states *within*
        # the target sector, and it can return two states that overlap
        # (or a duplicate) without saying so -- this is the one way that
        # shows up in the answer rather than in a per-state warning.
        warnings.warn(
            "get_dynamical_correlator(submode=\"SECTOR\"): the poles carry "
            "%.4f of the spectral weight of B|gs>, i.e. more than all of "
            "it -- the %d states of sector %s are not mutually orthogonal "
            "to the precision this sum assumes, and some weight is "
            "double-counted. Raise nsweeps/maxm, or lower nex."
            % (info["captured"], info["nex"],
               _sector_str_any(info["target_sector"])))
    if not quiet and info["captured"] < weight_tol:
        warnings.warn(
            "get_dynamical_correlator(submode=\"SECTOR\"): the %d states "
            "of sector %s account for %.1f%% of the spectral weight of "
            "B|gs> -- the rest sits at higher energies this truncated "
            "Lehmann sum does not describe. Increase nex, or use "
            "submode=\"KPM\"/\"TD\" if the response is a broad continuum."
            % (info["nex"], _sector_str_any(info["target_sector"]),
               100 * info["captured"]))
    es = np.array(es)
    out = np.zeros(len(es), dtype=np.complex128)
    for w, dE in zip(weights, poles):
        out = out + w * (1.0 / (es - dE + 1j * delta)
                         - 1.0 / (es - dE - 1j * delta))
    out = 1j * out / (2.0 * np.pi)
    info = dict(info, energies=poles, weights=weights, delta=delta)
    if return_poles:
        return es, out, info
    return es, out


# -- Layer 2: the physics-facing spectral functions ----------------------
#
# Both of these are thin: two (or three) sector-resolved correlators and
# an assembly. They exist because "give me A(w)" should not require the
# caller to know that the hole part is a separate sector solve, nor to get
# the chemical-potential bookkeeping right by hand -- and because the
# assembly is where the two conventions that are easy to get wrong live
# (which pole goes to negative frequency, and what mu is measured from).


def _lorentzian_sum(es, poles, weights, delta):
    """The same broadening kernel every submode here uses: dcex.py's
    1j*(advanced-retarded)/(2 pi), i.e. weight * delta/pi/((w-w_n)^2+delta^2)."""
    out = np.zeros(len(es), dtype=np.complex128)
    for w, e in zip(weights, poles):
        out = out + w * (1.0 / (es - e + 1j * delta)
                         - 1.0 / (es - e - 1j * delta))
    return 1j * out / (2.0 * np.pi)


def _fermion_operators(self, spin):
    """(C, Cdag) operator lists for the requested spin channel."""
    if spin is None:
        for c, cd in (("C", "Cdag"),):
            if hasattr(self, c) and hasattr(self, cd):
                return getattr(self, c), getattr(self, cd)
        raise ValueError(
            "get_spectral_function: this chain exposes no C/Cdag operators "
            "-- it is not a fermionic chain. Use get_spin_spectral_function "
            "for a spin chain.")
    if spin in ("up", "dn"):
        c, cd = "C" + spin, "Cdag" + spin
        if not hasattr(self, c):
            raise ValueError(
                "get_spectral_function(spin=%s): this chain has no %s "
                "operator -- spin= only applies to a spinful fermionic "
                "chain." % (repr(spin), c))
        return getattr(self, c), getattr(self, cd)
    raise ValueError("get_spectral_function: spin must be None, \"up\" or "
                     "\"dn\", got %s" % repr(spin))


def spectral_function(self, i=0, j=None, spin=None, nex=20, delta=2e-2,
                      es=np.linspace(-4.0, 4.0, 800), shift="mu",
                      return_poles=False, **kwargs):
    """The single-particle spectral function of a fermionic chain,

        A_ij(w) = sum_n <0|c_i|n^{N+1}><n^{N+1}|c_j^dag|0> L(w - w_n^+)
                + sum_m <0|c_j^dag|m^{N-1}><m^{N-1}|c_i|0> L(w - w_m^-)

    with the particle poles at w_n^+ = E_n^{N+1} - E_0^N - mu and the hole
    poles at w_m^- = -(E_m^{N-1} - E_0^N) - mu, both computed from the
    eigenstates of their own particle-number sector (see this module's
    docstring).

    The chemical potential mu = (E_0^{N+1} - E_0^{N-1})/2 is subtracted by
    default (`shift="mu"`), which is the convention the formula above is
    normally written in: particle weight then sits at w>0 and hole weight
    at w<0, separated by the gap. This is not cosmetic -- in the raw
    energies both families can land on the same side of zero, for any
    Hamiltonian carrying an explicit chemical-potential term or attractive
    interactions, and reading such a plot as "particle above, hole below"
    would be wrong. `shift=None` returns the unshifted axis; either way
    `mu` and the raw pole energies come back in the info dict.

    Returns `(es, A)`, or `(es, A, info)` with `return_poles=True`.
    """
    if j is None:
        j = i
    C, Cdag = _fermion_operators(self, spin)
    part_e, part_w, part_info = sector_poles(
        self, name=(C[i], Cdag[j]), nex=nex, **kwargs)
    hole_e, hole_w, hole_info = sector_poles(
        self, name=(Cdag[j], C[i]), nex=nex, **kwargs)
    # mu from the two sectors' own ground-state energies, not from the
    # lowest *weighted* pole: the lowest state of a sector can perfectly
    # well have zero matrix element with c|0>, and mu is a property of the
    # spectrum, not of the operator being probed.
    e_plus = min(part_info["energies_absolute"])
    e_minus = min(hole_info["energies_absolute"])
    mu = 0.5 * (e_plus - e_minus)
    off = mu if shift == "mu" else 0.0
    if shift not in ("mu", None):
        raise ValueError("get_spectral_function: shift must be \"mu\" or "
                         "None, got %s" % repr(shift))
    es = np.array(es)
    poles = np.concatenate([part_e - off, -hole_e - off])
    weights = np.concatenate([part_w, hole_w])
    out = _lorentzian_sum(es, poles, weights, delta)
    info = {"mu": mu, "shift": shift, "poles": poles, "weights": weights,
            "particle": {"energies": part_e, "weights": part_w,
                         "captured": part_info["captured"],
                         "sector": part_info["target_sector"]},
            "hole": {"energies": hole_e, "weights": hole_w,
                     "captured": hole_info["captured"],
                     "sector": hole_info["target_sector"]},
            "reference_sector": part_info["reference_sector"],
            "e0": part_info["e0"], "gap": e_plus + e_minus - 2 * part_info["e0"]}
    self._spectral_function_info = info
    if return_poles:
        return es, out, info
    return es, out


def spin_spectral_function(self, i=0, j=None, nex=20, delta=2e-2,
                           es=np.linspace(-0.5, 6.0, 800),
                           return_poles=False, **kwargs):
    """The dynamical spin structure factor of one site pair, resolved by
    Sz channel:

        S^{+-}_ij(w)  from the Sz+2 sector   (name=(S^-_i, S^+_j))
        S^{-+}_ij(w)  from the Sz-2 sector   (name=(S^+_i, S^-_j))
        S^{zz}_ij(w)  from the Sz sector itself

    Returns `(es, S)` with `S = S^zz + (S^{+-} + S^{-+})/2`, the same
    combination as S^xx+S^yy+S^zz, and -- with `return_poles=True` -- the
    three channels separately in the info dict. Keeping them separate is
    half the point of the method: which Sz channel a feature lives in is a
    selection rule no broadened single-curve method can show.

    Note S^zz's target sector *is* the reference sector, so its state n=0
    is the ground state itself and the channel carries an elastic w=0 pole
    of weight <0|Sz_i|0><0|Sz_j|0>. That is the correct Lehmann sum (the
    exact ED reference contains it too) rather than the connected
    structure factor; `connected=True` subtracts it.
    """
    if j is None:
        j = i
    connected = kwargs.pop("connected", False)
    Sp = operatornames.name2MO("Sp", self)
    Sm = operatornames.name2MO("Sm", self)
    Sz = operatornames.name2MO("Sz", self)
    channels, info = {}, {}
    for label, (A, B) in (("pm", (Sm[i], Sp[j])),
                          ("mp", (Sp[i], Sm[j])),
                          ("zz", (Sz[i], Sz[j]))):
        e, w, ch = sector_poles(self, name=(A, B), nex=nex, **kwargs)
        if label == "zz" and connected:
            keep = np.abs(e) > 1e-9  # drop the elastic ground-state pole
            e, w = e[keep], w[keep]
        channels[label] = {"energies": e, "weights": w,
                           "captured": ch["captured"],
                           "sector": ch["target_sector"]}
        info.setdefault("reference_sector", ch["reference_sector"])
    es = np.array(es)
    curves = {k: _lorentzian_sum(es, v["energies"], v["weights"], delta)
              for k, v in channels.items()}
    total = curves["zz"] + 0.5 * (curves["pm"] + curves["mp"])
    info.update({"channels": channels, "curves": curves,
                 "connected": connected})
    self._spin_spectral_function_info = info
    if return_poles:
        return es, total, info
    return es, total
