"""Lindblad_Chain: the vectorized Lindbladian of a spin or boson chain, its
spectrum and its steady state (Phase 2 of docs/clik_mps_lindblad_plan.md).

Anchors: the closed-form spectrum of a decaying spin, the dense
superoperator kron(A, B^T) built from the physical operator matrices (no
code shared with the doubled-chain MultiOperator but the matrices
themselves), and trace preservation <<1|L = 0.
"""

import numpy as np
import pytest

from dmrgpy import bosonchain, fermionchain, spinchain
from dmrgpy.lindblad import Lindblad_Chain, transpose
from dmrgpy.pyitensor.mpsalgebra import inner, randomMPS
from dmrgpy.pyitensor.vectorize import identity_mps


def _sp(sc, i):
    return sc.Sx[i] + 1j * sc.Sy[i]


def _sm(sc, i):
    return sc.Sx[i] - 1j * sc.Sy[i]


def test_decaying_spin_closed_form():
    """H = h Sz, L = sqrt(g) S-: lambda = 0, -g/2 +- i h, -g, <Sz> = -1/2."""
    h, g = 0.8, 0.3
    sc = spinchain.Spin_Chain([2], itensor_version="python")
    sc.set_hamiltonian(h * sc.Sz[0])
    lc = Lindblad_Chain(sc, [np.sqrt(g) * _sm(sc, 0)])
    exact = np.array([0, -g / 2 + 1j * h, -g / 2 - 1j * h, -g])
    for mode in ("ED", "DMRG"):
        lam, res, obs = lc.get_spectrum(mode=mode, nk=10)
        for x in exact:
            assert np.min(np.abs(lam[:4] - x)) < 1e-8
        assert abs(obs(sc.Sz[0], 0) + 0.5) < 1e-8


def _dense_superoperator(sc, H, jumps):
    ed = sc.get_ED_obj()

    def dense(op):
        return ed.get_operator(op).toarray()
    Js = [dense(J) for J in jumps]
    D = Js[0].shape[0]
    I = np.eye(D)
    Hd = dense(H)
    out = -1j * (np.kron(Hd, I) - np.kron(I, Hd.T))
    for Jd in Js:
        JdJ = Jd.conj().T @ Jd
        out += (np.kron(Jd, Jd.conj()) - 0.5 * np.kron(JdJ, I)
                - 0.5 * np.kron(I, JdJ.T))
    return out


def _spin_case():
    sc = spinchain.Spin_Chain([2] * 3, itensor_version="python")
    H = sum(sc.Sx[i] * sc.Sy[i + 1] - 0.3 * sc.Sy[i] * sc.Sz[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(2)) + 0.4 * sc.Sy[1]
    jumps = [0.5 * _sp(sc, 0) * sc.Sz[1], (0.3 + 0.2j) * sc.Sy[2],
             0.7 * _sm(sc, 1)]
    return sc, H, jumps


def _spin_one_case():
    sc = spinchain.Spin_Chain([3, 3], itensor_version="python")
    H = sc.Sx[0] * sc.Sy[1] + 0.3 * sc.Sz[0] * sc.Sz[0]
    jumps = [0.5 * _sp(sc, 0), 0.4 * sc.Sy[1] * sc.Sz[1]]
    return sc, H, jumps


def _boson_case():
    bc = bosonchain.Bosonic_Chain(3, maxnb=[3, 3, 3], itensor_version="python")
    H = sum(-1.0 * (bc.Adag[i] * bc.A[i + 1] + bc.Adag[i + 1] * bc.A[i])
            for i in range(2)) + 0.5 * sum(bc.N[i] * (bc.N[i] - 1)
                                           for i in range(3))
    jumps = [(bc.Adag[1] + bc.Adag[0]) * (bc.A[1] - bc.A[0]), 0.3 * bc.A[2],
             0.2j * bc.N[0]]
    return bc, H, jumps


@pytest.mark.parametrize("case", [_spin_case, _spin_one_case, _boson_case])
def test_vectorized_lindbladian_is_the_superoperator(case):
    """Tr(L^k), k = 1..5, of the doubled-chain operator against the dense
    kron superoperator: exact invariants, unlike eigenvalues, which eig
    resolves only to ~1e-8 at the degeneracies of the spin-1 case."""
    sc, H, jumps = case()
    sc.set_hamiltonian(H)
    lc = Lindblad_Chain(sc, jumps)
    A = lc.doubled.get_ED_obj().get_operator(lc.L).toarray()
    B = _dense_superoperator(sc, H, jumps)
    for k in range(1, 6):
        a = np.trace(np.linalg.matrix_power(A, k))
        b = np.trace(np.linalg.matrix_power(B, k))
        assert abs(a - b) <= 1e-12 * max(1.0, abs(b))


def test_trace_is_preserved():
    sc, H, jumps = _spin_case()
    sc.set_hamiltonian(H)
    lc = Lindblad_Chain(sc, jumps)
    S = lc.doubled._session
    one = identity_mps(S.sites)
    assert abs(inner(one, one) - 2 ** 3) < 1e-12
    rho = randomMPS(S.sites, 4, rng=np.random.default_rng(0))
    assert abs(inner(one, S._mpo(lc.L.to_terms()), rho)) < 1e-12


def test_driven_xxz_steady_state_against_ed():
    """Boundary drive and bulk dephasing on 4 sites (dissipative gap
    0.144). Measured at the defaults: steady state 1e-6 off, the gap mode
    2.2e-4, the magnetization profile 2.6e-6."""
    n = 4
    sc = spinchain.Spin_Chain([2] * n, itensor_version="python")
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
    assert abs(lc.steady_state_vev(sc.Sz[0]) - obs_ed(sc.Sz[0], 0)) < 1e-5


def test_transpose_of_names():
    sc = spinchain.Spin_Chain([2] * 2, itensor_version="python")
    ed = sc.get_ED_obj()
    for op in (sc.Sy[0] * sc.Sx[1], _sp(sc, 0) * sc.Sz[1], 1j * sc.Sy[1]):
        a = ed.get_operator(transpose(op)).toarray()
        b = ed.get_operator(op).toarray().T
        assert np.max(np.abs(a - b)) < 1e-14


def test_refuses_fermions():
    fc = fermionchain.Fermionic_Chain(2, itensor_version="python")
    with pytest.raises(NotImplementedError):
        Lindblad_Chain(fc, [fc.C[0]], H=fc.N[0])


def test_seeded_sector_bose_hubbard_against_ed():
    """The paper's model on 3 sites, 2 bosons. Number conservation gives
    one steady state per sector; seeding |psi><psi| and |psi><phi| keeps the
    space in N=2. Without the global subspace expansion of the first steps,
    two-site TDVP never moved these product starts (the next-nearest-
    neighbor terms of the doubled chain). Measured: occupations 1e-7, the
    slowest pair 1e-5 with residual 3e-6."""
    L, N = 3, 2
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N + 1] * L,
                                  itensor_version="python")
    h = sum(-(bc.Adag[j + 1] * bc.A[j] + bc.Adag[j] * bc.A[j + 1])
            for j in range(L - 1))
    h = h + 0.5 * sum(bc.N[j] * (bc.N[j] - 1) for j in range(L))
    bc.set_hamiltonian(h)
    bc.maxm = 64
    jumps = [(bc.Adag[j + 1] + bc.Adag[j]) * (bc.A[j + 1] - bc.A[j])
             for j in range(L - 1)]
    lc = Lindblad_Chain(bc, jumps)
    lam, res, obs = lc.get_spectrum(nk=20, dt=0.5, start=([1, 1, 0], [0, 1, 1]))
    n_dmrg = [obs(bc.N[j], 0).real for j in range(L)]
    bc.set_conserved_sector(Nb=N)
    lam_ed, _, obs_ed = lc.get_spectrum(mode="ED")
    assert abs(lam[0]) < 1e-6
    for j in range(L):
        assert abs(n_dmrg[j] - obs_ed(bc.N[j], 0).real) < 1e-5
    target = lam_ed[1] if lam_ed[1].imag > 0 else lam_ed[2]  # Im > 0 member
    k = int(np.argmin(np.abs(lam - target)))
    assert abs(lam[k] - target) < 1e-4
    assert res[k] < 1e-4


def _bose_hubbard(L, N, kappa=1.0):
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N + 1] * L,
                                  itensor_version="python")
    h = sum(-(bc.Adag[j + 1] * bc.A[j] + bc.Adag[j] * bc.A[j + 1])
            for j in range(L - 1))
    bc.set_hamiltonian(h + 0.5 * sum(bc.N[j] * (bc.N[j] - 1)
                                     for j in range(L)))
    bc.maxm = 64
    jumps = [np.sqrt(kappa) * (bc.Adag[j + 1] + bc.Adag[j])
             * (bc.A[j + 1] - bc.A[j]) for j in range(L - 1)]
    return bc, Lindblad_Chain(bc, jumps)


def test_sector_identity_is_the_left_steady_state_of_its_sector():
    """|1_q>> counts the configurations of charge q, and <<1_q|L = 0 when
    the jumps conserve the boson number; a decay breaks it."""
    from dmrgpy.pyitensor.chain import _is_left_zero_mode
    from dmrgpy.pyitensor.vectorize import sector_identity_mps
    bc, lc = _bose_hubbard(3, 2)
    S = lc.doubled._session
    one = sector_identity_mps(S.sites, 2)
    assert abs(inner(one, one) - 6) < 1e-12     # 6 ways to put 2 on 3
    assert abs(inner(one, identity_mps(S.sites)) - 6) < 1e-12
    assert _is_left_zero_mode(S._mpo(lc.L.to_terms()), one)
    lc2 = Lindblad_Chain(bc, [0.3 * bc.A[0]])
    assert not _is_left_zero_mode(S._mpo(lc2.L.to_terms()), one)


def test_stiff_lindbladian_takes_tdvp_substeps():
    """Ten times the dissipation of the seeded Bose-Hubbard test at
    dt=0.5. One TDVP step per snapshot is unstable there (backward bond
    evolution amplifies errors by exp(tau R)): measured 0.86 off ED in
    the occupations, the slowest pair 2 off. The trace-metered real
    branch doubles the sub-steps to 8: occupations 3e-4, the slowest pair
    1.2e-3 and 6.9e-3."""
    bc, lc = _bose_hubbard(3, 2, kappa=10.0)
    lam, res, obs = lc.get_spectrum(nk=20, dt=0.5, alpha=0.2,
                                    start=([1, 1, 0], [0, 1, 1]))
    n_dmrg = [obs(bc.N[j], 0).real for j in range(3)]
    bc.set_conserved_sector(Nb=2)
    lam_ed, _, obs_ed = lc.get_spectrum(mode="ED")
    assert abs(lam[0]) < 1e-6
    for j in range(3):
        assert abs(n_dmrg[j] - obs_ed(bc.N[j], 0).real) < 1e-3
    assert min(abs(lam - lam_ed[1])) < 2e-2
    assert min(abs(lam - lam_ed[2])) < 2e-2
