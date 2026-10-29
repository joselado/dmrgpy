"""get_krylov_spectrum: slowest eigenvalues of a generator from a Krylov
space of time-evolved MPS (pyitensor/clik.py, arXiv:2509.07709), Phase 1 of
docs/clik_mps_lindblad_plan.md.

Anchors: the dense spectrum of the same generator (mode="ED"), Cauchy
interlacing for a Hermitian H, and the residual of the reconstructed
eigenvector computed directly as an MPS. Tolerances are the measured
errors at the defaults (nk=20, dt=0.5, gram_tol=1e-12) with a margin of
about two.
"""

import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import applyMPO, inner, randomMPS, sum_many

from test_nh_dmrg import (nh_asymmetric_hopping_chain, nh_fermion_chain,
                          nh_pt_spin_chain)


def _heisenberg(n=6):
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    h = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    sc.set_hamiltonian(h)
    return sc, h


def test_hermitian_ritz_energies_interlace_and_converge():
    sc, _ = _heisenberg()
    ev, res, obs = sc.get_krylov_spectrum()
    eve, res_ed, obs_ed = sc.get_krylov_spectrum(mode="ED")
    E, Ee = (-ev).real, (-eve).real
    # Rayleigh-Ritz on a subspace: the k-th Ritz energy is never below the
    # k-th exact one, whatever the subspace
    assert np.all(E - Ee[:len(E)] >= -1e-9)
    assert abs(E[0] - Ee[0]) < 1e-9         # measured 1e-11
    assert abs(E[1] - Ee[1]) < 1e-7         # measured 1e-8
    assert res[0] < 2e-5                    # measured 5.5e-6
    assert np.max(res_ed[:3]) < 1e-12
    O = sc.Sz[0] * sc.Sz[1]
    assert abs(obs(O, 0) - obs_ed(O, 0)) < 1e-5   # measured 2.5e-7


def test_a_level_degenerate_beyond_the_two_starts_is_found_twice():
    """Two random starts span at most two copies of a degenerate level:
    the S=1 triplet above the Heisenberg ground state comes back twice."""
    sc, _ = _heisenberg()
    ev, _, _ = sc.get_krylov_spectrum()
    E = (-ev).real
    assert abs(E[1] - E[2]) < 1e-6
    assert E[3] > E[1] + 0.5


@pytest.mark.parametrize("model", [nh_asymmetric_hopping_chain,
                                   nh_pt_spin_chain, nh_fermion_chain])
def test_non_hermitian_lowest_levels_against_ed(model):
    sc, _ = model(6, "python")
    ev, _, _ = sc.get_krylov_spectrum()
    eve, _, _ = sc.get_krylov_spectrum(mode="ED")
    for lam in eve[:4]:
        assert np.min(np.abs(ev[:12] - lam)) < 5e-4   # measured <= 2.6e-4


def _true_residual(ks, G, k):
    c, lam = ks.coefficients[:, k], ks.eigenvalues[k]
    parts = []
    for ci, rho in zip(c, ks.states):
        g = applyMPO(G, rho, cutoff=0.0, maxdim=None)
        g.noPrime("Site")
        parts += [g * ci, rho * (-lam * ci)]
    r = sum_many(parts)
    return np.sqrt(abs(inner(r, r)))


def test_coefficient_residual_is_the_true_residual():
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    S = sc._session
    G = S._mpo((-1 * h).to_terms())
    ks = clik.krylov_spectrum(G, S.sites)
    for k in range(4):
        true = _true_residual(ks, G, k)
        assert abs(ks.residuals[k] - true) < 1e-3 * true + 1e-8


def test_mpo_product_order():
    n = 5
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
    A = sum(sc.Sx[i] * sc.Sx[i + 1] for i in range(n - 1)) + 0.7j * sc.Sz[0]
    B = sum(sc.Sz[i] * sc.Sy[i + 1] for i in range(n - 1)) + 0.3 * sc.Sx[2]
    S = sc._session
    Am, Bm = S._mpo(A.to_terms()), S._mpo(B.to_terms())
    rng = np.random.default_rng(1)
    x, y = randomMPS(S.sites, 4, rng=rng), randomMPS(S.sites, 4, rng=rng)
    By = applyMPO(Bm, y, cutoff=0.0, maxdim=None)
    By.noPrime("Site")
    assert abs(inner(x, clik.mpo_product(Am, Bm), y) - inner(x, Am, By)) < 1e-12


def test_refuses_what_it_cannot_do():
    sc, _ = _heisenberg(4)
    sc.set_conserved_sector(Sz=0)
    with pytest.raises(NotImplementedError):
        sc.get_krylov_spectrum()


def test_same_spectrum_on_the_jax_backend():
    pytest.importorskip("jax")
    from dmrgpy.pyitensor import backend as bk
    sc, h = nh_asymmetric_hopping_chain(6, "python")
    ev_np, _, _ = sc.get_krylov_spectrum(nk=8)
    try:
        bk.set_backend("jax")
        sc2, _ = nh_asymmetric_hopping_chain(6, "python")
        ev_jax, _, _ = sc2.get_krylov_spectrum(nk=8)
    finally:
        bk.set_backend("numpy")
    assert np.max(np.abs(ev_np[:4] - ev_jax[:4])) < 1e-8
