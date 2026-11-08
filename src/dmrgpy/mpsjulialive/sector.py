"""Conserved-sector mode on itensor_version="julia_live".

The same Many_Body_Chain.set_conserved_sector(Nf=6)/(Sz=0)/(Nf=8,Sz=0)
API as itensor_version=3 and "python", implemented with ITensors.jl's own
QN-carrying site indices (sites_qn.jl), so the confinement is structural
as on v3: ITensorMPS's dmrg() conserves the flux of its start state
exactly, and no penalty is needed.

The chain keeps two site sets. `jlsites_dense` is the one sites.jl builds,
the chain's reference, and never changes. `jlsites` is what every MPO and
every start state is built on: the QN set while a sector is set, the dense
set otherwise. So mpo.get_MPO, mps.start_mps and excited.jl build on the
right one without knowing about sectors, and promote_to_dense()/clearing
the sector only point `jlsites` back at the dense set.

The charge bookkeeping (which names a site type offers, which basis state
carries which charge, the in-sector product state of the start, the
normalization and charge check of every operator) is pyitensor's
(pyitensor/sector.py), run over a pyitensor SiteX built from the same site
codes and used for nothing else. The two agree on units and on basis-state
order (see sites_qn.jl), so the plan it returns is the start state here
too, and an operator that breaks the sector raises the same ValueError on
both backends, before ITensors ever sees it.
"""

from ..pyitensor import sector as pysector
from ..pyitensor.sites.siteset import SiteX, site_qn_names


def _qns(MBO):
    """The chain's sector as a sorted list of (name,target) pairs"""
    sector = getattr(MBO, "conserved_sector", None)
    return sorted((str(k), int(v)) for k, v in (sector or {}).items())


def _bookkeeping_sites(MBO, names):
    """pyitensor site set graded by `names`, for charges only. Raises the
    same ValueError pyitensor.chain.Chain.set_conserved_sector raises for a
    name no site carries, or a site that carries none of them."""
    codes = list(MBO.sites)
    for name in names:
        if not any(name in site_qn_names(code) for code in codes):
            raise ValueError(
                "set_conserved_sector: no site in this chain carries the "
                "quantum number %s (this chain offers %s)"
                % (name, ", ".join(sorted(set(
                    n for code in codes for n in site_qn_names(code))))
                   or "none at all"))
    return SiteX(codes, conserved=names)


def apply(MBO):
    """Put the chain's conserved sector on its Julia site set, or take it
    off. Returns True when a sector is now on the site set.

    Every check runs before anything changes, so a request that fails
    (an unknown name, an empty sector, a site type with no QN version on
    this backend) leaves the chain on the site set it was on."""
    qns = _qns(MBO)
    if not qns:
        MBO.jlsites = MBO.jlsites_dense
        MBO._jl_sector_plan = None
        return False
    if (getattr(MBO, "_jl_sector_plan", None) is not None
            and getattr(MBO, "_jl_sector_qns", None) == qns):
        # already on this sector's site set: keep it, and with it the Index
        # identity of every state built on it. A rejected request puts the
        # previous sector back through here, and a fresh QN set would leave
        # the chain's state on indices its operators no longer have.
        return True
    names = [k for k, _ in qns]
    if len(set(names)) != len(names):
        raise ValueError("set_conserved_sector: repeated quantum number "
                         "in %s" % (names,))
    sites = _bookkeeping_sites(MBO, names)
    plan = pysector.sector_state_plan(sites, qns)  # ValueError if empty
    for code in set(MBO.sites):
        if not (code in (0, 1, 2, 3) or 100 < code < 200):
            raise NotImplementedError(
                "set_conserved_sector on itensor_version=\"julia_live\" "
                "supports spin-1/2, spin-1, spinless and spinful fermion "
                "and boson sites; this chain has site code %s. Use "
                "itensor_version=3 or \"python\", or mode=\"ED\"." % code)
    from .juliasession import Main, to_julia_strvec
    from .sites import sites_text
    jl = Main.get_sites_qn(to_julia_strvec(sites_text(MBO)),
                           "Nf" in names, "Sz" in names, "Nb" in names)
    MBO.jlsites = jl
    MBO._jl_sector_plan = [int(s) for s in plan]
    MBO._jl_sector_qns = qns
    return True


def active(MBO):
    """True when the chain's Julia site set carries a sector"""
    return bool(_qns(MBO)) and getattr(MBO, "_jl_sector_plan", None) is not None


def start_state(MBO, linkdims):
    """Random MPS of link dimension `linkdims` inside the chain's sector
    (a Julia MPS), or None outside sector mode"""
    if not active(MBO):
        return None
    from .juliasession import Main
    return Main.sector_start_state(MBO.jlsites, start_plan(MBO), int(linkdims))


def start_plan(MBO):
    """The in-sector product state, as a Julia Vector{Int}, or None"""
    if not active(MBO):
        return None
    from juliacall import convert as jlconvert
    from .juliasession import Main
    return jlconvert(Main.Vector[Main.Int], list(MBO._jl_sector_plan))


def prepare_operator(MO, MBO):
    """The operator to build on the chain's site set: MO itself outside
    sector mode; inside one, MO rewritten with Sx/Sy expanded into S+/S-
    and exactly cancelling strings dropped, after checking that every
    remaining term conserves the sector (ValueError naming it otherwise,
    rather than a JuliaError from inside OpSum)."""
    if MBO is None or not active(MBO):
        return MO
    from ..multioperator import MultiOperator, MO2list
    qns = _qns(MBO)
    sites = _bookkeeping_sites(MBO, [k for k, _ in qns])
    terms = []
    for o in MO2list(MO):  # un-transformed names, 1-based sites
        factors = [(o[2 + 2 * k], int(o[3 + 2 * k]))
                   for k in range((len(o) - 2) // 2)]
        terms.append((complex(o[0], o[1]), factors))
    terms = pysector.sector_terms(sites, qns, terms)
    out = MultiOperator(term=False)
    for c, factors in terms:
        out.new_term(c=c)
        for name, site in factors:
            out.add_operator(name, site - 1)
    return out


def promote(MBO, jlmps):
    """A Julia MPS on the sector's QN sites, converted exactly to the
    chain's dense site set; returned unchanged when it is already dense"""
    from .juliasession import Main
    if not bool(Main.mps_hasqns(jlmps)):
        return jlmps
    return Main.promote_dense_mps(jlmps, MBO.jlsites_dense)
