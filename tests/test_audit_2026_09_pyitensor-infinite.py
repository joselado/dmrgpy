"""Regression tests for the 2026-09 audit's pyitensor infinite-chain cluster.

Findings #3, #8, #9, #20 and #32 of `docs/audit_2026_09_hole_hunt.md`. One
of them is a correctness bug and four are pure optimizations, so the tests
come in two shapes:

* #3 (`vumps_ms`'s sequential VUMPS returning energies BELOW the exact
  variational minimum, with `converged=True`) is checked against an exactly
  solvable model, and against the variational principle itself -- an energy
  below the exact minimum for a normalized state is impossible, so it is
  asserted as a one-sided bound as well as a two-sided tolerance.
* #8/#9/#20/#32 must not move a single number. Each is checked by running
  the optimized routine against the contraction it replaced, on the same
  inputs, and asserting agreement at ~1e-12 -- which is what an exact
  re-association of a contraction is allowed to differ by, and nothing
  more. Their *speed* is not asserted (a timing test on a shared box is a
  flake); the measured figures live in the audit entry and in the code
  comments at each fix.
"""

import numpy as np
import pytest

from dmrgpy import infinitechain
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor import idmrg_excitations as idmrg_exc
from dmrgpy.pyitensor import idmrg_window
from dmrgpy.pyitensor import vumps_ms


# ------------------------------------------------------------------ #3 --
# H = -FIELD * sum_i Sz_i + J * sum_i Sz_i Sz_{i+r} is diagonal in the Sz
# product basis (a classical Ising chain in a field), so its ground state is
# the fully polarized product state and its energy density is exact and
# closed-form at EVERY bond dimension >= 1. That is what makes it the right
# model here: any returned value below it is not "less converged", it is
# wrong -- and it is redundant bond dimension (the state needs D=1) that
# triggers the bug.
_FIELD, _J = 4.0, 0.7
_EXACT_POLARIZED = -_FIELD / 2.0 + _J / 4.0          # -1.825


def _polarized_chain(n_uc, reach, D, nrestarts=3):
    ic = infinitechain.Infinite_Spin_Chain(["1/2"] * n_uc)
    ic.maxm, ic.maxiter, ic.etol = D, 300, 1e-10
    ic.vumps_nrestarts = nrestarts
    h = 0
    for i in range(n_uc):
        h = h - _FIELD * ic.SzC[i]
        k = i + reach
        other = (ic.SzC[k] if k < n_uc
                 else ic.get_operator("Sz", k % n_uc, group=k // n_uc))
        h = h + _J * ic.SzC[i] * other
    ic.set_hamiltonian(h)
    return ic


@pytest.mark.parametrize("n_uc,reach,D", [
    (3, 1, 4),      # n_uc>2 alone routes to the sequential solver
    (3, 1, 6),
    (1, 2, 6),      # reach>1 routes there at any n_uc
    (4, 1, 6),
])
def test_sequential_vumps_never_returns_below_the_variational_minimum(
        n_uc, reach, D):
    """#3: `vumps_ms.ground_state` used to report `converged=True` at an
    energy far BELOW the exact minimum (-2.30, -1.8285, -4382 and -1.1e8
    against an exact -1.825), or to raise "every attempt at D=... failed",
    whenever `D` exceeded the bond dimension the state actually needs.

    The mechanism was `_cell_fixed_points`: a bare `eigs(k=1)` on a
    transfer matrix with a decoupled unimodular block returns an arbitrary
    element of a degenerate eigenspace, and normalizing a near-traceless
    one by its own trace blew the environment up. The state itself stayed
    exactly right (`<Sz>` reads 0.5 in every bad run), so only the energy
    read-off was wrong and nothing downstream flagged it.
    """
    ic = _polarized_chain(n_uc, reach, D)
    e0 = ic.gs_energy()
    # The one-sided bound is the real assertion: a normalized state cannot
    # have an energy below the exact minimum, at any bond dimension.
    assert e0 >= _EXACT_POLARIZED - 1e-6, (
        "energy {} is below the exact variational minimum {}".format(
            e0, _EXACT_POLARIZED))
    assert e0 == pytest.approx(_EXACT_POLARIZED, abs=1e-6)
    assert complex(ic.vev("Sz", 0)).real == pytest.approx(0.5, abs=1e-6)


def test_bond_fixed_points_are_exact_for_an_embedded_product_state():
    """#3, deterministically: the exactly-D=1 polarized state embedded in
    D=6 is the worst case -- its cell transfer matrix is the identity, so
    EVERY direction is a dominant eigenvector and the eigensolver has no
    well-posed answer at all. `_cell_fixed_points` must then return the
    fixed points the state itself names.

    Built by hand rather than by running a solve, since ARPACK's start
    vector is a randomness source the test cannot seed.
    """
    n_uc, D, d = 3, 6, 2
    # AL[l,s,r] = delta_{lr} delta_{s,0}: a legitimate isometry (every
    # site in the |up> state, bond basis carried straight through), whose
    # transfer is the identity map on (D,D) matrices.
    A = np.zeros((D, d, D), dtype=complex)
    for i in range(D):
        A[i, 0, i] = 1.0
    AL = [A.copy() for _ in range(n_uc)]
    AR = [A.copy() for _ in range(n_uc)]
    C = np.zeros((D, D), dtype=complex)
    C[0, 0] = 1.0                      # the state lives in one direction

    r_AL, l_AR = vumps_ms._cell_fixed_points(AL, AR, D, C_cell=C)
    expect = np.zeros((D, D), dtype=complex)
    expect[0, 0] = 1.0
    assert np.max(np.abs(r_AL - expect)) == pytest.approx(0.0, abs=1e-12)
    assert np.max(np.abs(l_AR - expect)) == pytest.approx(0.0, abs=1e-12)

    # Only the fixed points are asserted here, not a full `environments`
    # call on these tensors: this hand-built state is MORE degenerate than
    # anything a real solve produces (its cell transfer is exactly the
    # identity, so `I - T_cell + P` is rank 1 and the regularized
    # environment solve is genuinely singular). The end-to-end behaviour on
    # states a solver actually reaches is what
    # `test_sequential_vumps_never_returns_below_the_variational_minimum`
    # covers; this test pins the piece that can be checked without ARPACK's
    # unseedable start vector in the loop.


def test_zero_trace_fixed_point_is_rejected_rather_than_normalized():
    """#3, the normalization half: the old guard was `abs(tr) > 1e-300`,
    which let a trace of 1e-16 through and multiplied the environment by
    1e16. A candidate with no usable trace must raise instead."""
    M = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)   # traceless
    with pytest.raises(RuntimeError, match="zero"):
        vumps_ms._trace_normalized_hermitian(M, "test")


# ------------------------------------------------------------------ #9 --

def test_op_transfer_matrix_is_contiguous_and_unchanged():
    """#9: `_op_transfer_matrix` returned a stride-permuted (non-C-
    contiguous) view, so every consumer's `.reshape(chi*chi, -1)` silently
    copied the whole chi^4 tensor -- on EVERY application of a tensor that
    is built once and applied thousands of times. The values must be
    bit-for-bit what the transposed view held; only the layout changes.
    """
    rng = np.random.default_rng(0)
    D, d = 7, 3
    ket = rng.normal(size=(D, d, D)) + 1j * rng.normal(size=(D, d, D))
    bra = rng.normal(size=(D, d, D)) + 1j * rng.normal(size=(D, d, D))
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    for op in (None, M):
        E4 = idmrg_exc._op_transfer_matrix(ket, bra, op)
        assert E4.flags["C_CONTIGUOUS"], "the reshape in every consumer copies"
        k = ket if op is None else np.einsum('io,lir->lor', op, ket)
        reference = np.einsum('lpr,LpR->lLrR', k, np.conj(bra))
        assert np.max(np.abs(E4 - reference)) == pytest.approx(0.0, abs=1e-13)
        # and the reshape every consumer does is now a view, not a copy
        assert np.shares_memory(E4.reshape(D * D, -1), E4)


# ----------------------------------------------------------------- #20 --

def test_matrix_free_site_transfer_matches_the_rank4_application():
    """#20: `_dominant_fixed_point`'s ARPACK matvec applied materialized
    chi^4 transfer tensors although the same file's `_apply_site_transfer`
    does it in O(chi^3 d). The re-association is exact, so both directions
    must reproduce the rank-4 application to rounding -- including with a
    MIXED bra (what `imps_overlap` needs, and the reason the helpers grew
    an explicit `bra=` argument)."""
    rng = np.random.default_rng(1)
    D, d = 6, 3
    A = rng.normal(size=(D, d, D)) + 1j * rng.normal(size=(D, d, D))
    B = rng.normal(size=(D, d, D)) + 1j * rng.normal(size=(D, d, D))
    X = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    for bra in (None, B):
        E4 = idmrg_exc._op_transfer_matrix(A, A if bra is None else bra, None)
        assert np.max(np.abs(
            idmrg._apply_site_transfer(A, None, X, bra=bra)
            - idmrg._apply_transfer(E4, X))) == pytest.approx(0.0, abs=1e-12)
        assert np.max(np.abs(
            idmrg._apply_site_transfer_from_left(A, None, X, bra=bra)
            - idmrg._apply_transfer_from_left(E4, X))) == pytest.approx(
                0.0, abs=1e-12)


def test_fixed_points_are_identical_with_and_without_the_matrix_free_matvec():
    """#20, end to end: the `sites=` route through `_dominant_fixed_point`
    is a contraction-order change and nothing else, so both families of
    fixed points -- and the eigenvalue -- must come out the same."""
    ic = infinitechain.Infinite_Spin_Chain(["S=1/2"] * 2)
    ic.set_hamiltonian(ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1]
                       + ic.SzC[0] * ic.SzC[1]
                       + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                       + ic.SzC[1] * ic.SzR[0])
    ic.maxm, ic.maxiter = 16, 15
    ic.gs_method = "idmrg"
    ic.gs_energy()
    result = ic._result
    cell, n_cell = idmrg._correlator_cell(result)
    Es = idmrg._transfer_matrices(cell, n_cell)
    arrays = [idmrg._to_array_lpr(T) for T in cell]
    sites = (arrays, arrays)
    for fn in (idmrg._all_right_fixed_points, idmrg._all_left_fixed_points):
        plain = fn(Es, n_cell)
        free = fn(Es, n_cell, sites=sites)
        for a, b in zip(plain[0], free[0]):
            assert np.max(np.abs(a - b)) == pytest.approx(0.0, abs=1e-10)
        assert abs(plain[1] - free[1]) == pytest.approx(0.0, abs=1e-10)


# ------------------------------------------------------------- #8, #32 --

def _compose_and_close(bra_arrays, ket_arrays, l, rho_R):
    """`_close_array_chain`'s pre-fix inner contraction, verbatim: compose
    the whole chain into one rank-4 transfer tensor, then close it. Kept
    here as the reference the propagating form is checked against."""
    E = None
    for K, B in zip(ket_arrays, bra_arrays):
        step = np.einsum('lir,LiR->lLrR', K, np.conj(B))
        E = step if E is None else np.einsum('lLrR,rRsS->lLsS', E, step)
    return np.einsum('rR,rR->', np.einsum('lL,lLrR->rR', l, E), rho_R)


def test_propagated_chain_closure_matches_the_composed_one():
    """#8: `_close_array_chain` composed the full (chi,chi,chi,chi)
    transfer chain -- O(n chi^6) -- although the result was immediately
    closed on both ends, which propagating the boundary matrix does in
    O(n chi^3 d). Same scalar, only the contraction order differs.

    The ket and bra are given DIFFERENT bond dimensions on purpose: a
    shifted overlap between two independently evolved windows has exactly
    that shape, and it is what a naively-square implementation would break
    on."""
    rng = np.random.default_rng(2)
    d, nsites = 3, 7
    chi_k, chi_b = 5, 4
    ket = [rng.normal(size=(chi_k, d, chi_k))
           + 1j * rng.normal(size=(chi_k, d, chi_k)) for _ in range(nsites)]
    bra = [rng.normal(size=(chi_b, d, chi_b))
           + 1j * rng.normal(size=(chi_b, d, chi_b)) for _ in range(nsites)]
    l = rng.normal(size=(chi_k, chi_b)) + 1j * rng.normal(size=(chi_k, chi_b))
    rho = rng.normal(size=(chi_k, chi_b)) + 1j * rng.normal(size=(chi_k, chi_b))
    reference = _compose_and_close(bra, ket, l, rho)
    got = idmrg_window._propagate_close(bra, ket, l, rho)
    assert abs(got - reference) <= 1e-10 * abs(reference)


def test_window_environment_is_built_once_and_reused():
    """#32: `_close_array_chain` re-solved the call-invariant transfer-
    matrix fixed points on every invocation, and `local_expectation` built
    another copy of its own before calling it twice -- 25 full rebuilds for
    a 4-step `td_dynamical_correlator`. They are functions of the converged
    `IDMRGResult` alone, so they belong on it.

    Counted rather than timed: the number of `_transfer_matrices` calls a
    run makes is the thing the fix changes."""
    ic = infinitechain.Infinite_Spin_Chain(["S=1/2"] * 2)
    ic.set_hamiltonian(ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1]
                       + ic.SzC[0] * ic.SzC[1]
                       + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                       + ic.SzC[1] * ic.SzR[0])
    ic.maxm, ic.maxiter = 8, 12
    ic.gs_method = "idmrg"
    ic.gs_energy()

    calls = []
    original = idmrg._transfer_matrices

    def counted(*a, **k):
        calls.append(1)
        return original(*a, **k)

    idmrg_window._idmrg_mod._transfer_matrices = counted
    try:
        ks, ws, skw = ic.td_dynamical_correlator(
            "Sz", 0, "Sz", n_window=6, dt=0.2, nt=3, maxdim=8,
            x_values=[-1, 0, 1])
    finally:
        idmrg_window._idmrg_mod._transfer_matrices = original
    # One build for the window environment; the pre-fix code made one per
    # _close_array_chain call (18 here) plus local_expectation's own.
    assert len(calls) <= 2, "the window environment is being rebuilt per call"
    assert np.isfinite(skw).all() and np.max(np.abs(skw)) > 0.0


def test_window_local_expectation_is_uniform_on_the_ground_state():
    """#8/#32 together, on the public surface: with no perturbation applied
    the window's own `<Sz>` profile must read the ground state's (0 by
    symmetry on a Heisenberg chain) at every site. This is the check that
    would catch a re-association or a cache that changed the answer rather
    than only the cost."""
    ic = infinitechain.Infinite_Spin_Chain(["S=1/2"] * 2)
    ic.set_hamiltonian(ic.SxC[0] * ic.SxC[1] + ic.SyC[0] * ic.SyC[1]
                       + ic.SzC[0] * ic.SzC[1]
                       + ic.SxC[1] * ic.SxR[0] + ic.SyC[1] * ic.SyR[0]
                       + ic.SzC[1] * ic.SzR[0])
    ic.maxm, ic.maxiter = 8, 12
    ic.gs_method = "idmrg"
    ic.gs_energy()
    result = ic._result
    window = idmrg_window.build_window(result, 4)
    vals = [idmrg_window.local_expectation(window, result, s, "Sz")
            for s in range(1, window.mps.length() + 1)]
    assert np.max(np.abs(vals)) == pytest.approx(0.0, abs=1e-6)
