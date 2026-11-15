"""The levels of smallest real part of a non-Hermitian Hamiltonian by
imaginary time (nhlevels.py): Many_Body_Chain.get_nh_levels() and the
nh_gs_method="imaginary_time" route of gs_energy().

The PT-symmetric XX chain with i g (-1)^i Sz_i is free fermions after
Jordan-Wigner, so its levels are exact: the smallest real part is the sum
of the single-particle energies with negative real part, shared by every
subset of the purely imaginary ones. On 8 sites that is one pair, so four
levels share it: E0 twice (in two Sz sectors) and E0 +- i gamma. The other
models are tests/test_nh_dmrg.py's, against ED; the asymmetric hopping
chain is the one with H^T != H, which is what tests the left state."""

import numpy as np
import pytest

from dmrgpy import cppext, nhlevels, spinchain

from _helpers import setup_backend
from test_nh_dmrg import nh_asymmetric_hopping_chain, nh_fermion_chain

VERSIONS = [
    pytest.param(3, id="v3", marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    pytest.param("python", id="python"),
]

G = 0.3


def _pt_chain(n, version):
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    setup_backend(sc, version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
    for i in range(n):
        h = h + 1j * G * (-1)**i * sc.Sz[i]
    sc.set_hamiltonian(h)
    return sc, h


def _pt_levels(n):
    """The distinct levels sharing the smallest real part, and E0."""
    t = np.diag(1j * G * (-1.0)**np.arange(n)) \
        + 0.5 * (np.eye(n, k=1) + np.eye(n, k=-1))
    e = np.linalg.eigvals(t)
    e0 = e[e.real < -1e-12].sum().real
    gam = np.abs(e[np.abs(e.real) < 1e-9].imag).max()
    return np.array([e0, e0 + 1j * gam, e0 - 1j * gam]), e0


def _ed(sc):
    return np.linalg.eigvals(np.array(sc.get_ED_obj().get_hamiltonian().todense()))


@pytest.fixture
def no_fallback(monkeypatch):
    """Fail the test if the imaginary-time route hands over to NH-DMRG, so
    a pass means the route itself answered."""
    from dmrgpy import nhdmrg
    def refuse(*args, **kwargs):
        raise AssertionError("the imaginary-time route fell back to NH-DMRG")
    monkeypatch.setattr(nhdmrg, "gs_energy_nhdmrg", refuse)
    monkeypatch.setattr(nhdmrg, "nhdmrg", refuse)


@pytest.mark.parametrize("version", VERSIONS)
def test_get_nh_levels_returns_every_distinct_level_of_the_tied_real_part(version):
    """One NH-DMRG solve returns one of the four; the levels route returns
    the three distinct values, and no Ritz value below the spectrum passes
    its residual test."""
    np.random.seed(3)
    sc, _ = _pt_chain(8, version)
    exact, e0 = _pt_levels(8)
    es, rs = sc.get_nh_levels()
    low = es[es.real < e0 + 1e-3]
    assert len(low) >= 3
    for x in exact:
        assert np.min(np.abs(low - x)) < 1e-6
    assert np.all(es.real > e0 - 1e-6)
    assert np.all(rs < nhlevels.RITZ_TOL)
    assert np.all(np.diff(es.real) >= 0)


@pytest.mark.parametrize("version", VERSIONS)
def test_get_nh_levels_vectors_are_right_eigenvectors(version):
    np.random.seed(5)
    sc, h = _pt_chain(8, version)
    es, rs, vs = sc.get_nh_levels(vectors=True)
    assert len(vs) == len(es)
    for e, v in list(zip(es, vs))[:3]:
        r = h * v - e * v
        assert abs(r.dot(r))**0.5 / abs(v.dot(v))**0.5 < 1e-4


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("model", ["asymmetric", "fermion"])
def test_imaginary_time_ground_pair_matches_ed(version, model, no_fallback):
    """gs_energy() on the new route: the level of smallest real part, a
    right state and a left state that satisfy their own equations, the
    left one H^dagger|psil> = conj(E)|psil> on a non-symmetric H, paired
    with the stored right state."""
    build = {"asymmetric": nh_asymmetric_hopping_chain,
             "fermion": nh_fermion_chain}[model]
    np.random.seed(11)
    sc, h = build(6, version)
    ev = _ed(sc)
    sc.nh_gs_method = "imaginary_time"
    e = sc.gs_energy()
    assert e.real == pytest.approx(ev.real.min(), abs=1e-6)
    assert np.min(np.abs(ev - e)) < 1e-6
    psir, psil = sc.wf0, sc.nh_left_wf
    assert sc._nh_left_for is psir
    assert abs(psir.dot(psir) - 1) < 1e-10
    r = h * psir - e * psir
    l = h.get_dagger() * psil - np.conj(e) * psil
    assert abs(r.dot(r))**0.5 < 1e-4
    assert abs(l.dot(l))**0.5 / abs(psil.dot(psil))**0.5 < 1e-4
    assert sc._nh_last_residual[0] < 1e-4


@pytest.mark.parametrize("version", VERSIONS)
def test_imaginary_time_pair_on_the_tied_pt_chain(version, no_fallback):
    """The 10-site PT chain, four levels with one real part, where NH-DMRG
    used to stall: the stored pair is one of the three distinct levels."""
    np.random.seed(2)
    sc, h = _pt_chain(10, version)
    exact, e0 = _pt_levels(10)
    sc.nh_gs_method = "imaginary_time"
    e = sc.gs_energy()
    assert np.min(np.abs(exact - e)) < 1e-6
    assert sc._nh_last_residual[0] < 1e-4


@pytest.mark.parametrize("s,c", [(1e-6, 0.0), (1.0, 1e4)],
                         ids=["small-units", "offset"])
def test_imaginary_time_route_is_scale_covariant_and_offset_blind(s, c, no_fallback):
    """E(s*H + c) = s*E(H) + c on the new route: the grid is in units of
    the largest coefficient off the identity, so an offset does not shrink
    it and a small unit stretches it."""
    np.random.seed(13)
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    ev = _ed(sc)
    sc.set_hamiltonian(s * h + c)
    sc.nh_gs_method = "imaginary_time"
    e = sc.gs_energy()
    x = (e - c) / s
    assert np.min(np.abs(ev - x)) < 1e-6
    assert x.real == pytest.approx(ev.real.min(), abs=1e-6)
    assert sc._nh_last_residual[0] < 1e-4


def test_nh_gs_method_is_a_setting_and_part_of_the_solver_key(monkeypatch):
    """A constructor keyword like every setting, and a change of it between
    two bare gs_energy() calls solves again instead of returning the other
    route's stored answer."""
    sc = spinchain.Spin_Chain(["S=1/2"] * 4, itensor_version="python",
                              nh_gs_method="imaginary_time")
    assert sc.nh_gs_method == "imaginary_time"
    sc.nh_gs_method = "nhdmrg"
    h = 0
    for i in range(3):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
    sc.set_hamiltonian(h + 0.3j * sc.Sz[0])
    sc.gs_energy()
    calls = []
    real = nhlevels.ground_pair
    monkeypatch.setattr(nhlevels, "ground_pair",
                        lambda *a, **k: calls.append(1) or real(*a, **k))
    sc.nh_gs_method = "imaginary_time"
    sc.gs_energy()
    assert calls == [1]


def test_unknown_nh_gs_method_raises():
    sc, _ = _pt_chain(4, "python")
    sc.nh_gs_method = "imaginary"
    with pytest.raises(ValueError, match="nh_gs_method"):
        sc.gs_energy()


@pytest.mark.skipif(not cppext.available(2), reason="requires the v2 extension")
def test_backends_without_the_route_use_nhdmrg():
    """v2 has no implicit steps: the setting falls through to NH-DMRG."""
    sc, _ = _pt_chain(6, 2)
    assert not nhlevels.available(sc)
    ev = _ed(sc)
    sc.nh_gs_method = "imaginary_time"
    e = sc.gs_energy()
    assert e.real == pytest.approx(ev.real.min(), abs=1e-6)


def _count_fallbacks(monkeypatch):
    from dmrgpy import nhdmrg
    calls = []
    real = nhdmrg.gs_energy_nhdmrg
    monkeypatch.setattr(nhdmrg, "gs_energy_nhdmrg",
                        lambda *a, **k: calls.append(1) or real(*a, **k))
    return calls


def test_the_better_pair_is_kept_when_both_miss(monkeypatch, capsys):
    """When NH-DMRG does no better than a pair above the tolerance, the
    imaginary-time pair is kept, marked unconverged, with a warning."""
    from dmrgpy import nhdmrg
    real_pair = nhlevels.ground_pair
    def fake_nhdmrg(self, **kwargs):
        e, l, r = real_pair(self, self.hamiltonian)
        self._nh_last_residual = (1.0, 1e-4)
        return e + 0.5, l, r
    monkeypatch.setattr(nhdmrg, "nhdmrg", fake_nhdmrg)
    np.random.seed(29)
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    ev = _ed(sc)
    e = sc.gs_energy(tol=1e-30)  # nothing passes
    assert np.min(np.abs(ev - e)) < 1e-6
    assert sc._gs_unconverged is not None
    assert "kept over NH-DMRG" in capsys.readouterr().out


def test_a_route_that_raises_falls_back_to_nhdmrg(monkeypatch):
    """A failure inside the route, a LinAlgError from the small
    eigenproblems included, hands the solve to NH-DMRG."""
    calls = _count_fallbacks(monkeypatch)
    def broken(*args, **kwargs):
        raise np.linalg.LinAlgError("Eigenvalues did not converge")
    monkeypatch.setattr(nhlevels, "ground_pair", broken)
    np.random.seed(17)
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    ev = _ed(sc)
    e = sc.gs_energy()
    assert calls == [1]
    assert e.real == pytest.approx(ev.real.min(), abs=1e-6)


def test_a_pair_that_fails_the_certificate_falls_back_to_nhdmrg(monkeypatch):
    """A pair that is not an eigenpair is not stored: NH-DMRG runs, and
    its better pair is the one kept."""
    from dmrgpy import nhdmrg
    calls = []
    real = nhdmrg.nhdmrg
    monkeypatch.setattr(nhdmrg, "nhdmrg",
                        lambda *a, **k: calls.append(1) or real(*a, **k))
    def wrong(self, H=None, **kwargs):
        r = self.random_state()
        return -1.0 + 0j, r.copy(), r
    monkeypatch.setattr(nhlevels, "ground_pair", wrong)
    np.random.seed(19)
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    ev = _ed(sc)
    e = sc.gs_energy()
    assert calls == [1]
    assert e.real == pytest.approx(ev.real.min(), abs=1e-6)


def test_get_nh_levels_on_a_chain_pinned_to_ed():
    """The route evolves an MPS even when the chain's own random state is
    an ED vector."""
    np.random.seed(23)
    sc, _ = _pt_chain(6, "python")
    sc.mode = "ED"
    exact, e0 = _pt_levels(6)
    es, rs = sc.get_nh_levels()
    assert np.min(np.abs(es - exact[0])) < 1e-6
