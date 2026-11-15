"""Regressions of docs/audit_2026_11_15_hole_hunt.md, finding 3: nhdmrg()
on v2, v3 and "python" settled on a sweep fixed point that is not an
eigenpair when four eigenvalues share the smallest real part, already on
the 10-site PT-symmetric XX chain (relative residual 0.29 to 0.44, Re E up
to 0.03 below the spectrum, on every attempt).

The cure is in the sweep (pyitensor/nhdmrg.py's module docstring): the
left state is carried in the reference's transposed convention, so the
shared isometry of a complex-symmetric H holds r's own Schmidt basis, and
the local solver restarts on the 3/5 most wanted Ritz vectors
(Krylov-Schur) instead of the single selected one. The reference value is
exact: the chain is free fermions after Jordan-Wigner (hopping 1/2, onsite
i g (-1)^i, and a constant -i g/2 sum_i (-1)^i), so its smallest real part
is the sum of the single-particle eigenvalues with negative real part.

Each test runs one attempt (ntries=1), so the retry loop cannot hide a
stalled draw, from three seeds of the global numpy stream; v2 and v3 draw
their start inside ITensor, so the seeds vary only the "python" start and
the three v2/v3 runs are three independent draws."""

import numpy as np
import pytest

from dmrgpy import cppext, spinchain
from dmrgpy.pyitensor import nhdmrg as pnh

from _helpers import setup_backend

VERSIONS = [
    pytest.param(2, marks=pytest.mark.skipif(
        not cppext.available(2),
        reason="requires the compiled mpscpp2 (ITensor v2) extension")),
    pytest.param(3, marks=pytest.mark.skipif(
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
    sc.maxm, sc.nsweeps, sc.cutoff = 30, 10, 1e-10
    return sc, h


def _free_fermion_levels(n):
    """(smallest real part of the many-body spectrum, single-particle
    eigenvalues) of the PT chain, from its Jordan-Wigner free fermions."""
    t = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        t[i, i + 1] = t[i + 1, i] = 0.5
    for i in range(n):
        t[i, i] = 1j * G * (-1)**i
    e = np.linalg.eigvals(t)
    const = -0.5j * G * sum((-1)**i for i in range(n))
    return (e[e.real < -1e-12].sum() + const).real, e


def test_free_fermion_reference_matches_ed():
    """The reference itself, against full ED of the 8-site chain, and the
    degeneracy that makes the chain hard: one pair of modes with zero real
    part, so four many-body levels share the smallest real part."""
    n = 8
    sc, _ = _pt_chain(n, "python")
    ev = np.linalg.eigvals(np.array(sc.get_ED_obj().get_hamiltonian().todense()))
    emin, e1 = _free_fermion_levels(n)
    assert ev.real.min() == pytest.approx(emin, abs=1e-10)
    assert np.sum(np.abs(ev.real - emin) < 1e-8) == 4
    assert np.sum(np.abs(e1.real) < 1e-9) == 2


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("seed", [1234, 1, 2])
def test_nhdmrg_reaches_an_eigenpair_on_the_ten_site_pt_chain(version, seed):
    n = 10
    emin, e1 = _free_fermion_levels(n)
    # the four levels of smallest real part: emin, emin, emin +- i gamma
    gamma = np.abs(e1[np.abs(e1.real) < 1e-9].imag).max()
    levels = np.array([emin, emin + 1j * gamma, emin - 1j * gamma])
    np.random.seed(seed)
    sc, h = _pt_chain(n, version)
    e, psil, psir = sc.nhdmrg(ntries=1)
    assert np.min(np.abs(levels - e)) < 1e-6
    assert e.real > emin - 1e-6  # the failure sat up to 0.03 below
    assert psil.dot(psir) == pytest.approx(1.0, abs=1e-10)
    r = h * psir - e * psir
    l = h.get_dagger() * psil - np.conj(e) * psil
    assert abs(r.dot(r))**0.5 / abs(psir.dot(psir))**0.5 < 1e-4
    assert abs(l.dot(l))**0.5 / abs(psil.dot(psil))**0.5 < 1e-4


SELS = [("SR", 0.0), ("closest", 0.3 - 0.2j), ("SRTieBreak", -1.0 + 0.5j)]


@pytest.mark.parametrize("sel,target", SELS)
def test_ritz_order_leads_with_select_ritz(sel, target):
    """The thick restart keeps a leading slice of _ritz_order, which has to
    start with the pair _select_ritz returns, ties included."""
    rng = np.random.default_rng(7)
    for trial in range(200):
        m = int(rng.integers(1, 12))
        ev = rng.normal(size=m) + 1j * rng.normal(size=m)
        if m > 2 and trial % 3 == 0:  # a real-part tie, as on the PT chain
            ev[1] = ev[0].real + 1j * rng.normal()
            ev[2] = ev[0].real - 1j * ev[1].imag
        order = pnh._ritz_order(ev, sel, target)
        assert sorted(order.tolist()) == list(range(m))
        assert order[0] == pnh._select_ritz(ev, sel, target)


def test_krylov_schur_resolves_a_cluster_of_equal_real_parts():
    """A non-normal 300x300 operator whose four eigenvalues of smallest real
    part share one real part, two of them a conjugate pair, as on the PT
    chain: from a random start the thick restart returns an eigenpair of
    that real part, at krylovdim=20 with 6 restarts."""
    rng = np.random.default_rng(3)
    d = 300
    lam = np.concatenate([[-2.0, -2.0, -2.0 + 0.3j, -2.0 - 0.3j],
                          rng.uniform(-1.5, 3.0, d - 4)
                          + 1j * rng.uniform(-1.0, 1.0, d - 4)])
    S = np.eye(d) + 0.2 * rng.normal(size=(d, d)) / np.sqrt(d)
    A = S @ np.diag(lam) @ np.linalg.inv(S)
    x0 = rng.normal(size=d) + 1j * rng.normal(size=d)
    e, x = pnh._krylov_schur_smallest_real(lambda v: A @ v, x0, 20, 6)
    x = x / np.linalg.norm(x)
    assert e.real == pytest.approx(-2.0, abs=1e-8)
    assert np.linalg.norm(A @ x - e * x) < 1e-7
