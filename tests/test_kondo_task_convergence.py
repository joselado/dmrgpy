"""Which poles of the TASK Kondo spectrum the Krylov space has converged
(pyitensor/task.py's residual, kondospectrumtk/taskkondo.py's convergence).

At T=0 every pole is a sharp step of dI/dV. A Ritz state with a residual
is a quadrature node, not an eigenstate, and its energy moves with nkry and
between two solves of the same chain (3e-3 of the peak of dI/dV on a
24-site chain at nkry=100, at the biases a pole had crossed). The residual
is the measure of it, and get_kondo_spectrum warns.
"""

import warnings

import numpy as np
import pytest

from dmrgpy import cppext, taskdmrg
from dmrgpy.kondospectrumtk.taskkondo import TaskKondoSpectrum
from dmrgpy.pyitensor import task

from test_kondo_spectrum_task import EV, KW, impurity_chain
from test_task_batched_spectrum import _tensors

VERSIONS = ["python", pytest.param(3, marks=pytest.mark.skipif(
    not cppext.available(3), reason="v3 extension not compiled"))]


@pytest.mark.parametrize("real", [False, True])
@pytest.mark.parametrize("nkry", [2, 5])
def test_the_residual_is_the_norm_of_what_h_takes_out_of_the_basis(real, nkry):
    rand, Ms, Ws = _tensors(7, real=real)
    ts = task.TangentSpace(Ms, Ws)
    starts = [ts.gauge(ts.unflat(rand(ts.flat_size()))) for _ in range(3)]
    store, n, Hq, (j0, R) = task.block_krylov(ts, starts, nkry, residual=True)
    ws, U = np.linalg.eigh(Hq)
    Ul = U[j0:, :]
    res = np.sqrt(np.abs(np.einsum("jf,jk,kf->f", Ul.conj(), R, Ul)))
    for f in (0, n // 2, n - 1):
        v = task.ritz_vector((ts, store, n), U[:, f])
        hv = ts.apply_h(v)
        explicit = ts.norm([a - ws[f] * b for a, b in zip(hv, v)])
        assert res[f] == pytest.approx(explicit, rel=1e-8, abs=1e-12 * abs(ws).max())
    assert np.max(res) > 1e-3 * np.max(np.abs(ws))  # a few blocks converge nothing


def test_an_exhausted_krylov_space_has_no_residual():
    rand, Ms, Ws = _tensors(8, d=2)
    Ms, Ws = Ms[:4], Ws[:4]
    Ms[-1] = Ms[-1][:, :, :1] if Ms[-1].shape[2] > 1 else Ms[-1]
    Ms[3] = rand(Ms[2].shape[2], 2, 1)
    Ws[3] = Ws[3][:, :, :, :1]
    ts = task.TangentSpace(Ms, Ws)
    starts = [ts.gauge(ts.unflat(rand(ts.flat_size())))]
    _, n, _, (j0, R) = task.block_krylov(ts, starts, 10 * ts.dimension(),
                                         residual=True)
    assert n == ts.dimension() and j0 == n and R.size == 0


@pytest.mark.parametrize("version", VERSIONS)
def test_a_saturated_spectrum_is_converged_and_does_not_warn(version):
    sc = impurity_chain(version=version)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **KW)
    c = TaskKondoSpectrum(sc, 0, nkry=50).convergence(np.max(np.abs(EV)))
    assert c["poles"] > 0 and c["unconverged"] == 0
    assert c["weight"] == 0.0 and c["step"] == 0.0


@pytest.mark.parametrize("version", VERSIONS)
def test_a_short_krylov_space_warns_and_says_how_much_is_open(version):
    sc = impurity_chain(version=version)
    with pytest.warns(RuntimeWarning, match="are not converged at nkry=3"):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=3, **KW)
    ks = TaskKondoSpectrum(sc, 0, nkry=3)
    c = ks.convergence(np.max(np.abs(EV)))
    assert 0 < c["unconverged"] <= c["poles"]
    assert 0.0 < c["weight"] <= 1.0 and 0.0 < c["step"] <= 1.0
    assert ks.residual[0] == 0.0 and len(ks.residual) == ks.dim
    # nothing above the window is counted
    assert ks.convergence(0.0)["poles"] == 0


def test_every_pole_is_within_its_residual_of_an_eigenvalue():
    """What the residual bounds: the distance from a Ritz energy to the
    nearest eigenvalue of P H P, here those of the exhausted Krylov space.
    A pole with no residual is an eigenvalue, the others are not yet."""
    sc = impurity_chain()
    small = TaskKondoSpectrum(sc, 0, nkry=8)
    exact = TaskKondoSpectrum(sc, 0, nkry=50)
    assert np.max(exact.residual) < 1e-9 * np.max(np.abs(exact.e))
    width = np.max(np.abs(exact.e))
    dist = np.array([np.min(np.abs(exact.e - e)) for e in small.e])
    assert np.all(dist <= small.residual + 1e-9 * width)
    moved = dist > 1e-6 * width
    assert moved.any() and np.all(small.residual[moved] > 1e-6 * width)


def test_the_residual_is_optional_in_the_session_spectrum():
    sc = impurity_chain()
    taskdmrg.ground_state(sc)
    ops = [(sc.Sz[0].to_terms(), sc.Sz[0].get_dagger().to_terms())]
    assert len(taskdmrg.session_spectrum(sc, ops, 5)) == 3
    e, mats, dim, res = taskdmrg.session_spectrum(sc, ops, 5, residual=True)
    assert res.shape == e.shape and res[0] == 0.0 and np.all(res >= 0.0)
    res[:] = 7.0  # a copy, as the rest
    assert not np.any(taskdmrg.session_spectrum(sc, ops, 5, residual=True)[3] == 7.0)
