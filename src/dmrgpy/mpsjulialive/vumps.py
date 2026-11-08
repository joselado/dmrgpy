"""Bridge from Infinite_Many_Body_Chain(..., itensor_version="julia_live")
to the native Julia VUMPS solver in vumps.jl (module DmrgpyVumps).

What runs where. The Hamiltonian automaton, a dense array W[a,i,o,b] for
the grouped unit cell, is assembled on the Python side by the same
pyitensor/idmrg.py builder the "python" backend uses (`_build_automaton`,
then `vumps._group_automaton`), and so are the dense operator matrices the
observables need (site matrices, the Jordan-Wigner string, the kron
embedding into the grouped site). This is data preparation only: the VUMPS
iteration (fixed points, environments, the H_AC/H_C eigenproblems, the
polar updates, the D-ramp with restarts) and every contraction with the
state run in Julia, and the state itself stays on the Julia side.

Scope: gs_method="vumps", n_uc in {1, 2}, couplings reaching at most one
unit cell (the grouped reach-1 path of pyitensor/vumps.py). Anything else
raises NotImplementedError in infinitechain.py before reaching here.
"""

import numpy as np

SUPPORTED = ('gs_method="vumps" with n_uc<=2 and nearest-cell couplings '
             '(gs_energy, vev, correlator)')


class JuliaVumpsResult:
    """The Julia VumpsState plus the Python-side site data the observables
    need to build their operator matrices."""

    def __init__(self, jl_state, sites_uc, n_uc):
        self.jl = jl_state
        self.sites_uc = sites_uc
        self.n_uc = n_uc
        self.e0 = float(jl_state.e0)
        self.converged = bool(jl_state.converged)
        self.gauge_mismatch = float(jl_state.mismatch)
        self.niter = int(jl_state.niter)


def _main():
    from .juliasession import Main
    return Main


def grouped_automaton(h_intra_op, h_inter_op, site_types, n_uc):
    """(sites_uc, W, hunit, reach_one) for the grouped cell, built exactly
    as pyitensor.vumps.vumps_ground_state builds them."""
    from ..pyitensor import idmrg
    from ..pyitensor import vumps as pyvumps
    from ..pyitensor import idmrg_excitations as idmrg_exc
    sites_uc, W_bulk = idmrg._build_automaton(h_intra_op, h_inter_op,
                                              site_types, n_uc)
    hunit = pyvumps._hamiltonian_unit(h_intra_op, h_inter_op, sites_uc, n_uc)
    W = pyvumps._group_automaton(W_bulk, n_uc)
    return sites_uc, np.asarray(W, dtype=complex), hunit, idmrg_exc.is_reach_one(W)


def ground_state(site_types, h_intra_op, h_inter_op, n_uc, D, tol, maxiter,
                 niter, nrestarts, verbose=False):
    if n_uc > 2:
        raise NotImplementedError(
            "Infinite_Many_Body_Chain on itensor_version=\"julia_live\": "
            "n_uc={} is not implemented, only {}".format(n_uc, SUPPORTED))
    if D < 1:
        raise ValueError("vumps: D must be >= 1, got {}".format(D))
    if nrestarts < 1:
        raise ValueError("vumps: nrestarts must be >= 1, got {}".format(nrestarts))
    sites_uc, W, hunit, reach_one = grouped_automaton(
        h_intra_op, h_inter_op, site_types, n_uc)
    if not reach_one:
        raise NotImplementedError(
            "Infinite_Many_Body_Chain on itensor_version=\"julia_live\": a "
            "coupling reaching further than one unit cell needs the "
            "sequential multi-site VUMPS, which is not implemented there, "
            "only {}".format(SUPPORTED))
    st = _main().DmrgpyVumps.vumps_ground_state(
        W, int(n_uc), int(D), float(tol), int(maxiter), int(niter),
        int(nrestarts), float(hunit), bool(verbose))
    return JuliaVumpsResult(st, sites_uc, n_uc)


def onsite_expectation(result, opname, p):
    """<opname> at sub-site p of the cell, read off AC in Julia."""
    from ..pyitensor import vumps as pyvumps
    M = pyvumps._embed_group_operator(
        result.sites_uc, result.n_uc,
        {p: result.sites_uc.site_type(p + 1).matrix(opname)})
    return complex(_main().DmrgpyVumps.onsite_expectation(
        result.jl, np.asarray(M, dtype=complex)))


def two_point_correlator(result, opname_i, p_i, opname_j, r):
    """<opname_i(p_i) opname_j(p_i+r)>, the operator matrices (with their
    Jordan-Wigner strings) built exactly as pyitensor.vumps.two_point_correlator
    builds them, the contraction with AC and AR done in Julia."""
    from ..pyitensor import idmrg
    from ..pyitensor import vumps as pyvumps
    if r < 0:
        raise ValueError("two_point_correlator: r must be >= 0")
    sites_uc, n_uc = result.sites_uc, result.n_uc
    embed = lambda ops: np.asarray(
        pyvumps._embed_group_operator(sites_uc, n_uc, ops), dtype=complex)
    cell_offset, p_j = divmod(p_i + r, n_uc)
    term = [1.0, [opname_i, p_i], [opname_j, p_i + r]]
    _rel, coef, mats, ferm = idmrg._term_site_matrices(term, sites_uc, n_uc)
    strung = ferm[0] if r > 0 else False

    def string_ops(lo, hi):
        if not strung:
            return {}
        return {p: sites_uc.site_type(p + 1).matrix("F") for p in range(lo, hi)}

    jl = _main().DmrgpyVumps
    if cell_offset == 0:
        if p_j == p_i:
            M = embed({p_i: mats[0]})
        else:
            ops = {p_i: mats[0], p_j: mats[1]}
            ops.update(string_ops(p_i + 1, p_j))
            M = embed(ops)
        return complex(coef * complex(jl.onsite_expectation(result.jl, M)))
    ops_i = {p_i: mats[0]}
    ops_i.update(string_ops(p_i + 1, n_uc))
    ops_j = {p_j: mats[1]}
    ops_j.update(string_ops(0, p_j))
    M_cross = embed(string_ops(0, n_uc)) if strung else None
    val = jl.two_point_correlator(result.jl, embed(ops_i), M_cross,
                                  embed(ops_j), int(cell_offset))
    return complex(coef * complex(val))
