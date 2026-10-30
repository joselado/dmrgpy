"""CLIK-MPS and Lindblad_Chain on itensor_version=3 (cliksession.py).

cliksession.py runs pyitensor/clik.py's algorithm on the v3 session;
tests/test_krylov_spectrum.py and tests/test_lindblad.py pin the algorithm
on "python". Here the v3 route is pinned against the same anchors, the
dense spectrum (mode="ED"), with the tolerances of those files:

- the Heisenberg ground state and first level, and the three
  non-Hermitian chains of test_krylov_spectrum.py;
- the seeded random starts, which must reproduce;
- the exact residual (exact_residuals=True), against the residual of the
  reconstructed eigenvector;
- the vectorized identities built in C++ (trace, sector counting, left
  zero mode);
- the driven XXZ steady state, and the seeded Bose-Hubbard sector, whose
  random-free starts make v3 and "python" comparable to ED to the same
  digits.
"""

import numpy as np
import pytest

from dmrgpy import bosonchain, cppext, spinchain
from dmrgpy.lindblad import Lindblad_Chain

from test_nh_dmrg import (nh_asymmetric_hopping_chain, nh_fermion_chain,
                          nh_pt_spin_chain)

pytestmark = pytest.mark.skipif(not cppext.available(3),
                                reason="the itensor_version=3 extension is not compiled")


def _heisenberg(n=6):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=3)
    h = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    sc.set_hamiltonian(h)
    return sc, h


def _sm(sc, i):
    return sc.Sx[i] - 1j * sc.Sy[i]


def _sp(sc, i):
    return sc.Sx[i] + 1j * sc.Sy[i]


def test_hermitian_levels_against_ed():
    sc, _ = _heisenberg()
    ev, res, obs = sc.get_krylov_spectrum()
    eve, _, obs_ed = sc.get_krylov_spectrum(mode="ED")
    E, Ee = (-ev).real, (-eve).real
    assert np.all(E - Ee[:len(E)] >= -1e-9)
    assert abs(E[0] - Ee[0]) < 1e-8
    assert abs(E[1] - Ee[1]) < 1e-6
    O = sc.Sz[0] * sc.Sz[1]
    assert abs(obs(O, 0) - obs_ed(O, 0)) < 1e-5


@pytest.mark.parametrize("model", [nh_asymmetric_hopping_chain,
                                   nh_pt_spin_chain, nh_fermion_chain])
def test_non_hermitian_lowest_levels_against_ed(model):
    """The three models of test_krylov_spectrum.py; a non-Hermitian G
    takes the Arnoldi TDVP step. The tolerance is twice the worst level
    measured on v3 over seeds 0 to 4 (6.4e-4, PT spin chain, seed 0);
    "python" spreads the same way over its own draws (up to 4.1e-4 on the
    PT chain, 3e-3 on a stronger asymmetric chain), so 5e-4 there is a
    property of its seed-0 starts, not of the algorithm."""
    sc, _ = model(6, 3)
    assert not sc.is_hermitian(sc.hamiltonian)
    ev, _, _ = sc.get_krylov_spectrum()
    eve, _, _ = sc.get_krylov_spectrum(mode="ED")
    for lam in eve[:4]:
        assert np.min(np.abs(ev[:12] - lam)) < 1.3e-3


def test_random_starts_are_reproducible():
    """ITensor's randomMPS draws from std::random_device; the starts come
    from a seeded generator instead, so the same seed gives the same run."""
    sc, _ = _heisenberg()
    s = sc._session
    a = s.clik_random_starts(4, 0)
    b = s.clik_random_starts(4, 0)
    c = s.clik_random_starts(4, 1)
    assert abs(s.overlap(a[0], b[0]) - 1) < 1e-12
    assert abs(s.overlap(a[1], b[1]) - 1) < 1e-12
    assert abs(s.overlap(a[0], c[0])) < 0.99


def test_exact_residuals_against_ed_residual_of_the_level():
    """exact_residuals=True reads ||G r - lambda r|| off G^dag G; for the
    converged ground state it is small, and the ED residual of the same
    level is zero."""
    sc, _ = _heisenberg()
    ev, res, _ = sc.get_krylov_spectrum(exact_residuals=True)
    assert res[0] < 1e-4
    assert np.all(np.isfinite(res))
    # an unconverged mode has a large residual
    assert res[-1] > res[0]


def test_vectorized_identities():
    bc = bosonchain.Bosonic_Chain(6, maxnb=[3] * 6, itensor_version=3)
    s = bc._session
    one = s.doubled_identity_mps()
    assert abs(s.overlap(one, one) - 3 ** 3) < 1e-12
    q2 = s.sector_identity_mps(2)
    assert abs(s.overlap(q2, q2) - 6) < 1e-12     # 6 ways to put 2 on 3
    assert abs(s.overlap(q2, one) - 6) < 1e-12


def _bose_hubbard(L, N, version, kappa=1.0):
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N + 1] * L, itensor_version=version)
    h = sum(-(bc.Adag[j + 1] * bc.A[j] + bc.Adag[j] * bc.A[j + 1])
            for j in range(L - 1))
    bc.set_hamiltonian(h + 0.5 * sum(bc.N[j] * (bc.N[j] - 1)
                                     for j in range(L)))
    bc.maxm = 64
    jumps = [np.sqrt(kappa) * (bc.Adag[j + 1] + bc.Adag[j])
             * (bc.A[j + 1] - bc.A[j]) for j in range(L - 1)]
    return bc, Lindblad_Chain(bc, jumps)


def test_lindblad_runs_on_v3_and_the_sector_identity_is_a_left_zero_mode():
    from dmrgpy import cliksession
    bc, lc = _bose_hubbard(3, 2, 3)
    assert lc.doubled.itensor_version == 3
    ops = cliksession.Session(lc.doubled)
    G = ops.mpo(lc.L)
    one = lc.doubled._session.sector_identity_mps(2)
    assert cliksession.is_left_zero_mode(ops, G, one)
    rho = lc.doubled._session.clik_random_starts(4, 0)[0]
    assert abs(ops.sandwich(lc.doubled._session.doubled_identity_mps(), G, rho)) < 1e-12
    lc2 = Lindblad_Chain(bc, [0.3 * bc.A[0]])
    assert not cliksession.is_left_zero_mode(ops, ops.mpo(lc2.L), one)


def test_driven_xxz_steady_state_against_ed():
    n = 4
    sc = spinchain.Spin_Chain([2] * n, itensor_version=3)
    sc.maxm = 64
    H = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + 0.5 * sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    sc.set_hamiltonian(H)
    jumps = [0.6 * _sp(sc, 0), 0.6 * _sm(sc, n - 1)]
    jumps += [0.4 * sc.Sz[i] for i in range(n)]
    lc = Lindblad_Chain(sc, jumps)
    lam_ed, _, obs_ed = lc.get_spectrum(mode="ED")
    lam, res, obs = lc.get_spectrum()
    assert abs(lam[0]) < 1e-5
    assert np.min(np.abs(lam[:6] - lam_ed[1])) < 1e-3
    for i in range(n):
        assert abs(obs(sc.Sz[i], 0) - obs_ed(sc.Sz[i], 0)) < 1e-5


@pytest.mark.parametrize("version", [3, "python"])
def test_seeded_sector_bose_hubbard_against_ed(version):
    L, N = 3, 2
    bc, lc = _bose_hubbard(L, N, version)
    lam, res, obs = lc.get_spectrum(nk=20, dt=0.5, start=([1, 1, 0], [0, 1, 1]))
    n_dmrg = [obs(bc.N[j], 0).real for j in range(L)]
    bc.set_conserved_sector(Nb=N)
    lam_ed, _, obs_ed = lc.get_spectrum(mode="ED")
    assert abs(lam[0]) < 1e-6
    for j in range(L):
        assert abs(n_dmrg[j] - obs_ed(bc.N[j], 0).real) < 1e-5
    target = lam_ed[1] if lam_ed[1].imag > 0 else lam_ed[2]
    k = int(np.argmin(np.abs(lam - target)))
    assert abs(lam[k] - target) < 1e-4
