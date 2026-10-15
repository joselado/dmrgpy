"""get_fidelity(fmode="derivative") on a degenerate ground state.

With a ground state degenerate within delta, the estimator returned
|det Ud|, the determinant of the second difference of the overlap matrix,
which scales as chi^g for g degenerate states. It now takes the second
difference of |det U|, the matrix form of the second difference of the
overlap used for a single state, which for g states is the trace of the
quantum metric of the manifold.
"""
import contextlib
import io

import numpy as np
import pytest

from dmrgpy import fidelity, spinchain
from dmrgpy.fidelity import _chi_manifold, _chi_single


def _quiet(f, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return f(*args, **kwargs)


def _ising(n):
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    h1 = 2*sum(sc.Sx[i] for i in range(n))
    sc.set_hamiltonian(h0 + h1)
    return sc, h0, h1


def _dense(sc, op):
    m = sc.get_ED_obj().get_operator(op)
    return m.toarray() if hasattr(m, "toarray") else np.asarray(m)


def test_unique_ground_state_scalar_and_matrix_forms_agree():
    n, lam, dl = 8, 1.2, 0.02
    sc, h0, h1 = _ising(n)
    small = _quiet(fidelity.get_fidelity, sc, h0, h1, lam, dl=dl,
                   delta=1e-4, mode="ED")
    large = _quiet(fidelity.get_fidelity, sc, h0, h1, lam, dl=dl,
                   delta=1e-1, mode="ED")
    # the same three ground states through the scalar and the 1x1 matrix
    states = []
    for x in (lam, lam + dl, lam - dl):
        sc.set_hamiltonian(h0 + x*h1)
        states.append(_quiet(sc.get_gs, mode="ED"))
    w, wp, wm = states
    single = _chi_single(w, wp, wm, dl)
    matrix = _chi_manifold([w], [wp], [wm], dl)
    # by hand, from full ED eigenvectors
    v = [np.linalg.eigh(_dense(sc, h0 + x*h1))[1][:, 0]
         for x in (lam, lam + dl, lam - dl)]
    hand = (2 - abs(v[0].conj() @ v[1]) - abs(v[0].conj() @ v[2]))/dl**2
    assert single == pytest.approx(matrix, rel=1e-9)
    assert small == pytest.approx(hand, rel=1e-6)
    assert large == pytest.approx(hand, rel=1e-6)
    assert single == pytest.approx(hand, rel=1e-6)


def test_degenerate_pair_gives_subspace_susceptibility():
    # lambda=0: the two Ising ground states are exactly degenerate
    n, dl = 8, 0.02
    sc, h0, h1 = _ising(n)
    chi = _quiet(fidelity.get_fidelity, sc, h0, h1, 0.0, dl=dl,
                 delta=1e-6, mode="ED")
    H0 = _dense(sc, h0)
    H1 = _dense(sc, h1)
    es, vs = np.linalg.eigh(H0)
    gs = np.abs(es - es[0]) < 1e-9
    assert np.sum(gs) == 2
    w = np.abs(vs[:, ~gs].conj().T @ H1 @ vs[:, gs])**2
    hand = np.sum(w/((es[~gs] - es[0])**2)[:, None])
    # one flip at an edge costs 2, in the bulk 4, for each of the pair
    assert hand == pytest.approx(1 + (n-2)/8, rel=1e-12)
    # finite-difference error of order dl^2, measured 2.5e-5 relative
    assert chi == pytest.approx(hand, rel=1e-3)
    # the old |det Ud| gave 0.766 here
    assert chi > 1.5
