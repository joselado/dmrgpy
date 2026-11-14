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


# ---- T>0 on the TASK route: the electrodes at T, the sample in its
# ground state, refused when that is not the thermal state

KB = 8.617333262e-5


def _cold(sc, factor=45.):
    """A temperature at which the lowest excitation has a Boltzmann factor
    of exp(-factor), so that mode="ED" occupies the ground state only."""
    e = np.sort(np.linalg.eigvalsh(sc.get_ED_obj().get_operator(
        sc.hamiltonian).todense()))
    gap = e[e > e[0] + 1e-9 * (e[-1] - e[0])][0] - e[0]
    return gap / (factor * KB)


@pytest.mark.parametrize("version", VERSIONS)
def test_thermal_electrodes_match_ed_while_the_sample_is_cold(version):
    sc = impurity_chain(version=version)
    kw = dict(KW, T=_cold(sc))
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **kw)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **kw)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))
    # and it is not the T=0 curve: the steps are broadened
    _, cold = sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **KW)
    assert np.max(np.abs(tk - cold)) > 1e-3 * np.max(np.abs(ed))


def test_thermal_electrodes_average_a_degenerate_manifold_like_ed():
    from dmrgpy import spinchain
    sc = spinchain.Spin_Chain(["S=1/2"] * 3, itensor_version="python")
    h = 0
    for i in range(2):
        h = h + 1e-3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + sc.Sz[i] * sc.Sz[i + 1])
    h = h + 0.3e-3 * sc.Sz[0] * sc.Sz[2]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 8, 20, 0.
    eV = np.linspace(-4e-3, 4e-3, 81)
    kw = dict(site=0, Jrho_s=-0.1, U=0.2, T=_cold(sc), omega0=20e-3,
              Gamma0=5e-5, order=3)
    _, ed = sc.get_kondo_spectrum(eV, mode="ED", **kw)
    _, tk = sc.get_kondo_spectrum(eV, mode="DMRG", submode="TASK", nkry=30,
                                  n_gs=2, **kw)
    assert np.max(np.abs(tk - ed)) < 1e-9 * np.max(np.abs(ed))


@pytest.mark.parametrize("version", VERSIONS)
def test_a_warm_sample_is_refused(version):
    """With the lowest excitation occupied, mode="ED" sums over it and the
    DMRG route cannot: it raises and names the route that can."""
    sc = impurity_chain(version=version)
    warm = dict(KW, T=_cold(sc, factor=2.))
    with pytest.raises(ValueError, match='use mode="ED"'):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="TASK", nkry=50, **warm)
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **warm)  # which does run
    assert np.all(np.isfinite(ed))


def test_the_bound_is_the_boltzmann_factor_of_the_lowest_pole():
    from dmrgpy.kondospectrumtk import taskkondo
    sc = impurity_chain()
    e1 = np.min(TaskKondoSpectrum(sc, 0, nkry=50).e[1:])
    edge = e1 / (KB * np.log(1. / taskkondo.THERMAL_OCCUPATION))
    TaskKondoSpectrum(sc, 0, nkry=50, T=0.99 * edge)
    with pytest.raises(ValueError, match="Boltzmann factor"):
        TaskKondoSpectrum(sc, 0, nkry=50, T=1.01 * edge)


def test_other_submodes_and_a_negative_temperature_are_refused():
    sc = impurity_chain()
    with pytest.raises(ValueError, match="only supports T=0"):
        sc.get_kondo_spectrum(EV, mode="DMRG", submode="KPM", delta=1e-5,
                              es=np.linspace(-1e-2, 1e-2, 50),
                              **dict(KW, T=0.1, order=2))
    with pytest.raises(ValueError, match="T must be >= 0"):
        sc.get_kondo_spectrum(EV, mode="DMRG", **dict(KW, T=-1.))


@pytest.mark.skipif(not cppext.available(3), reason="v3 extension not compiled")
def test_a_chain_too_short_for_the_session_is_the_thermal_ed_spectrum():
    """Two sites on v3 fall back to ED (mode.py), and there T is the
    temperature of everything, the sample included."""
    from dmrgpy import spinchain
    sc = spinchain.Spin_Chain(["S=1", "S=1/2"], itensor_version=3)
    h = 1e-3 * sc.Sz[0] * sc.Sz[0] + 2e-4 * sc.Sz[0] + 5e-4 * (
        sc.Sx[0] * sc.Sx[1] + sc.Sy[0] * sc.Sy[1] + sc.Sz[0] * sc.Sz[1])
    sc.set_hamiltonian(h)
    kw = dict(KW, T=2.0)  # warm: excited states are occupied
    _, ed = sc.get_kondo_spectrum(EV, mode="ED", **kw)
    _, tk = sc.get_kondo_spectrum(EV, mode="DMRG", **kw)
    assert np.max(np.abs(tk - ed)) < 1e-12 * np.max(np.abs(ed))


def test_a_thermal_step_wider_than_the_pole_spacing_needs_no_warning():
    sc = impurity_chain()
    ks = TaskKondoSpectrum(sc, 0, nkry=3)
    emax = np.max(np.abs(EV))
    open_ = ks.convergence(emax)["unconverged"]
    assert open_ > 0
    ks.T = 1e-9  # a step far narrower than any spacing: nothing changes
    assert ks.convergence(emax)["unconverged"] == open_
    ks.T = np.max(ks.e) / KB  # a step as wide as the whole spectrum
    assert ks.convergence(emax)["unconverged"] == 0
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ks.warn_unconverged(emax)


def _odd_chain(version, n=9, maxm=10):
    """An odd Heisenberg chain: a Kramers doublet, here at a bond dimension
    that truncates it, so the partner of a member is not at zero exactly."""
    from dmrgpy import spinchain
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + 1e-3 * (sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                        + sc.Sz[i] * sc.Sz[i + 1])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, 20
    return sc


@pytest.mark.parametrize("version", VERSIONS)
def test_a_truncated_manifold_is_still_averaged_at_finite_temperature(version):
    sc = _odd_chain(version)
    eV = np.linspace(-4e-3, 4e-3, 81)
    kw = dict(site=0, Jrho_s=-0.1, U=0.2, T=_cold(sc, factor=12.),
              omega0=20e-3, Gamma0=5e-5, order=3)
    _, ed = sc.get_kondo_spectrum(eV, mode="ED", **kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the split of the truncated doublet
        _, tk = sc.get_kondo_spectrum(eV, mode="DMRG", nkry=60, n_gs=2, **kw)
        # one state of a doublet is not the thermal state, and it says so
        with pytest.raises(ValueError, match="Pass n_gs=2"):
            sc.get_kondo_spectrum(eV, mode="DMRG", nkry=60, **kw)
    assert np.max(np.abs(tk - ed)) < 1e-3 * np.max(np.abs(ed))


def test_a_lowest_pole_that_is_not_converged_cannot_pass_the_bound():
    """A Ritz energy is an upper bound on the excitation, so a short
    Krylov space would let a warm sample through."""
    sc = impurity_chain()
    exact = np.min(TaskKondoSpectrum(sc, 0, nkry=50).e[1:])
    short = TaskKondoSpectrum(sc, 0, nkry=1)
    k = 1 + np.argmin(short.e[1:])
    assert short.e[k] > 1.05 * exact and short.residual[k] > 0.05 * exact
    # warm for the true excitation, cold for the pole of the short space
    T = 0.5 * (exact + short.e[k]) / (KB * np.log(1e3))
    with pytest.raises(ValueError, match="Boltzmann factor"):
        TaskKondoSpectrum(sc, 0, nkry=50, T=T)
    with pytest.raises(ValueError, match="is not converged at nkry=1"):
        TaskKondoSpectrum(sc, 0, nkry=1, T=T)


def test_a_state_below_the_one_handed_over_is_refused_at_finite_temperature(
        monkeypatch):
    from test_audit_2026_10_kondo import heisenberg
    e = np.array([0., -0.5, 0.3])
    M = np.diag([0.5, 0.5, -0.5]).astype(complex)
    M[0, 2] = M[2, 0] = 0.4
    M[1, 2] = M[2, 1] = 0.1  # the spin reaches the lower state
    monkeypatch.setattr(taskdmrg, "session_spectrum",
                        lambda chain, ops, nkry, residual=False, **kw: (
                            e, [M, 0 * M, M], 10, np.zeros(3), 0.0))
    sc = heisenberg(3, itensor_version="python")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        TaskKondoSpectrum(sc, 0, nkry=5)  # T=0 warns, as before
        with pytest.raises(ValueError, match="not the ground state"):
            TaskKondoSpectrum(sc, 0, nkry=5, T=1e-3, kB=1.)
